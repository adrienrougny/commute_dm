"""`momapy_bel` elements as real activity-flow CellDesigner content.

`commute_dm.bel_export` loads a BEL knowledge graph's influence-graph projection
as a `momapy_bel.core.BELModel`. This module is the half that knows what a BEL
*element* is worth in CellDesigner: it describes each element in CellDesigner
vocabulary and builds the frozen momapy species and glyphs. `bel_export` imports
`bel2cd`, never the reverse.

Everything is read off the elements' fields -- `modifications`, `variants`,
`fragment`, `location`, `members`, `abundance` -- so nothing here parses a BEL
string or touches the database.

What an element becomes::

    ProteinAbundance("GDE1")                      GenericProtein("GDE1")
    + ProteinModification("pho", "Y", "18")       + Modification(residue Y18, PHOSPHORYLATED)
    + Variant("P, T, 212")                        + Modification(residue T212, PHOSPHORYLATED)
    + Variant("K,670,N")                          + StructuralState("K670N")
    + Fragment("672_713")                         TruncatedProtein + StructuralState
    + Location(GO, "intracellular")               compartment=Compartment("intracellular")
    Activity(abundance=ProteinAbundance("CSNK2A1"))  GenericProtein("CSNK2A1", active=True)
    ComplexAbundance(members={A, B})              Complex("A:B", subunits={A, B})
    CompositeAbundance(members={A, B})            Complex + StructuralState("composite")
    Abundance(CHEBI, "dopamine")                  SimpleMolecule("dopamine")
    Abundance(CONSO, "Tau aggregates")            Unknown("Tau aggregates")
    BiologicalProcess / Pathology                 Phenotype

**Four identity invariants.** Naming a species by its whole BEL string used to
guarantee uniqueness by *value*; bare HGNC symbols do not, so uniqueness is
carried by *object identity* instead. Each is enforced by a check in
`commute_dm.submaps` that runs before every write, because getting one wrong
produces a file that writes without error and then fails to read back with a
`KeyError`.

(a) **One template object per `(template_class, name)` within the export.**
`submaps.make_submap_from_model_elements` derives `model.species_templates` *by
value*; two value-equal but distinct template objects collapse there while each
species keeps its own object, so the writer emits one `<protein>` and the other
species emit a `<proteinReference>` to an id that was never written. See
:func:`make_templates`. Identity with the *stored* collections' templates is
not established here: the export is standalone, and the save merges value-equal
elements into one database node (hash integration with a seeded
`object_key_to_node`).

(b) **`Modification.residue` *is* an object inside its species' template's
residue container.** The writer emits the modification's residue id and the
template's residue ids independently; after `renumber_ids` two value-equal
residue objects have two different ids and the modification names a residue that
is not in the list, which the reader looks up unguarded. Hence a template's
residues are the union over *every* proteoform of that protein, and the residue
lookup is re-derived from the *interned* template.

(c) **A subunit object is never the same object as a top-level species.** The
writer keys `build_subunit_to_complex` by `id(subunit)` and skips any species
whose `id()` is in that index. `ProteinAbundance("APP")` is both a projected
element and a member of several complexes, so memoising species by element --
the obvious way to get value-collapse -- would make APP's own glyph disappear.
:func:`make_species` therefore always builds a **fresh object tree**; interning
happens only at the top level, in `bel_export`.

This one is an **in-memory property of the export only, and does not survive the
save**: a top-level `p(X)` and the same `p(X)` inside a complex are value-equal,
so hash integration stores them as one node and hydration returns one object.
A source map built from the stored collection therefore has hundreds of subunits
that *are* their top-level species (651 in AD), and that is expected --
`submaps.make_submap_from_model_elements` promotes such a species to its own
object per sub-map, before the invariant checks run. Never run
`check_identity_invariants` on a source map; it is a pre-write check for
sub-maps.

(d) **One compartment object per location, hanging off an undrawn root.** Same
collapse as (a), one level down -- and the root parent is not decorative:
`pd2af.utils._build_dot_graph` only attaches a compartment's dot cluster to the
graph when `compartment.outside is not None` (there is no `else`), so a drawn
compartment with no `outside` gets an orphan cluster, its species never reach
graphviz, and they keep their throwaway positions. The root also keeps BEL
compartments structurally distinct from the stored CellDesigner ones, 18 of whose
names a BEL `loc()` shares: a BEL `nucleus` has `outside=<BEL root>`, a stored
one `outside=default`, so they can never be value-equal and neither can the
species inside them.
"""

import dataclasses
import math
import re

import momapy.celldesigner
import momapy.coloring
import momapy.core.layout
import momapy.drawing
import momapy.geometry
import momapy_bel.core
import pd2af.celldesigner.building_layout

import commute_dm.submaps


# ---------------------------------------------------------------------------
# describing an element in CellDesigner vocabulary
# ---------------------------------------------------------------------------

# The BEL element classes that have a fixed CellDesigner counterpart, mapped to
# `(species class, template class)`. `Activity` is absent -- it has no class of
# its own and takes its subject's -- and so is `Abundance`, which splits on its
# namespace (see below).
#
# This table is **not** cosmetic: it selects the template class too, so an entry
# decides both the glyph and whether the species declares a `<protein>` /
# `<gene>` / `<rna>`.
#
# A microRNA is RNA; `AntisenseRNA` would assert a strand the KG does not give.
# CellDesigner has no "composite", and a set of abundances is closest to a
# complex; it has no disease class either, and a pathology is a state of the
# system, so both `bp()` and `path()` are phenotypes.
BEL_CLASS_TO_CD_CLASSES = {
    momapy_bel.core.ProteinAbundance: (
        momapy.celldesigner.GenericProtein,
        momapy.celldesigner.GenericProteinTemplate,
    ),
    momapy_bel.core.GeneAbundance: (
        momapy.celldesigner.Gene,
        momapy.celldesigner.GeneTemplate,
    ),
    momapy_bel.core.RNAAbundance: (
        momapy.celldesigner.RNA,
        momapy.celldesigner.RNATemplate,
    ),
    momapy_bel.core.MicroRNAAbundance: (
        momapy.celldesigner.RNA,
        momapy.celldesigner.RNATemplate,
    ),
    momapy_bel.core.ComplexAbundance: (momapy.celldesigner.Complex, None),
    momapy_bel.core.CompositeAbundance: (momapy.celldesigner.Complex, None),
    momapy_bel.core.BiologicalProcess: (momapy.celldesigner.Phenotype, None),
    momapy_bel.core.Pathology: (momapy.celldesigner.Phenotype, None),
}

# Only a CHEBI abundance is reliably a small molecule (350 of the AD KG's 427:
# drugs and metabolites). The rest are not molecules at all and `SimpleMolecule`
# would misstate them -- MESH mixes proteins (`Acetylcholinesterase`), cell types
# (`Astrocytes`) and structures (`Blood-Brain Barrier`); CONSO holds aggregates
# (`Tau aggregates`); GO holds cellular components (`mitochondrion`). `Unknown`
# is CellDesigner's entity-of-unspecified-type, which is what they are to us.
ABUNDANCE_MOLECULE_NAMESPACES = frozenset({"CHEBI"})

# A fragment turns a protein into CellDesigner's truncated one. A gene that is
# also fragmented therefore yields two `<protein>` declarations (GENERIC and
# TRUNCATED), which is idiomatic CellDesigner.
_TRUNCATED_CLASSES = (
    momapy.celldesigner.TruncatedProtein,
    momapy.celldesigner.TruncatedProteinTemplate,
)

# Species classes by what they can *carry*, which is what the modifier routing
# below branches on.
PROTEIN_SPECIES_CLASSES = (
    momapy.celldesigner.GenericProtein,
    momapy.celldesigner.TruncatedProtein,
)
COMPLEX_SPECIES_CLASSES = (momapy.celldesigner.Complex,)

# Where a template keeps the residues a `Modification` can name, and which class
# those residues are. A `GeneTemplate` / `RNATemplate` keeps them in `regions` as
# `ModificationSite`s instead, which the writer emits differently -- gene and RNA
# badges are deferred, so enabling them later is one entry here.
TEMPLATE_RESIDUE_FIELDS = {
    momapy.celldesigner.GenericProteinTemplate: (
        "modification_residues",
        momapy.celldesigner.ModificationResidue,
    ),
    momapy.celldesigner.TruncatedProteinTemplate: (
        "modification_residues",
        momapy.celldesigner.ModificationResidue,
    ),
    momapy.celldesigner.GeneTemplate: (
        "regions",
        momapy.celldesigner.ModificationSite,
    ),
    momapy.celldesigner.RNATemplate: (
        "regions",
        momapy.celldesigner.ModificationSite,
    ),
}

# A BEL `pmod` type code to the CellDesigner state its badge shows. A
# `ProteinModification` with a namespace is an ontology pmod
# (`pmod(GO:"protein oxidation")`), which names a process rather than a residue
# state and becomes a structural state.
PMOD_TYPE_TO_MODIFICATION_STATE = {
    "pho": momapy.celldesigner.ModificationState.PHOSPHORYLATED,
    "ubi": momapy.celldesigner.ModificationState.UBIQUITINATED,
    "ace": momapy.celldesigner.ModificationState.ACETYLATED,
    "me0": momapy.celldesigner.ModificationState.METHYLATED,
    "gly": momapy.celldesigner.ModificationState.GLYCOSYLATED,
    "ogl": momapy.celldesigner.ModificationState.GLYCOSYLATED,
}

# A BEL1 modification code, as written inside a `var()`, to the CellDesigner
# state its badge shows -- `var(P,S,396)` is `pmod(Ph,Ser,396)` in BEL1 spelling,
# and 90 of the four knowledge graphs' 213 variants are one of these. The codes
# are the spec's `bel1_migration.protein_modifications` table; the ones with no
# CellDesigner state (`F` farnesylation, `R` ADP-ribosylation, `S` SUMOylation)
# are deliberately absent, and a code outside this table is not read as a
# modification at all -- which is also what keeps `S` free to be serine.
VARIANT_CODE_TO_MODIFICATION_STATE = {
    "P": momapy.celldesigner.ModificationState.PHOSPHORYLATED,
    "Ph": momapy.celldesigner.ModificationState.PHOSPHORYLATED,
    "A": momapy.celldesigner.ModificationState.ACETYLATED,
    "Ac": momapy.celldesigner.ModificationState.ACETYLATED,
    "M": momapy.celldesigner.ModificationState.METHYLATED,
    "Me": momapy.celldesigner.ModificationState.METHYLATED,
    "U": momapy.celldesigner.ModificationState.UBIQUITINATED,
    "Ub": momapy.celldesigner.ModificationState.UBIQUITINATED,
    "G": momapy.celldesigner.ModificationState.GLYCOSYLATED,
    "Glyco": momapy.celldesigner.ModificationState.GLYCOSYLATED,
    "H": momapy.celldesigner.ModificationState.HYDROXYLATED,
    "Hy": momapy.celldesigner.ModificationState.HYDROXYLATED,
}

# BEL `var()` descriptors are free text in practice. Two shapes are recognised,
# in this order, after whitespace is stripped (so `"P, S, 9"` and `"P,S,9"` are
# the same variant): a modification written as a variant, and an amino-acid
# substitution. **They are told apart by which field holds the number**: a
# modification is `(code, amino acid, position)` and a substitution is
# `(amino acid, position, amino acid)`, so `P,S,396` is phosphorylation of S396
# while `P,396,S` would be Pro396 -> Ser. The letters alone cannot decide it --
# `P` is both proline and phosphorylation. Anything else is kept verbatim as a
# structural state -- `misfolded`, `p.D614G`, `c.863G>A`, `del`.
_VARIANT_MODIFICATION = re.compile(r"^([A-Za-z]+),([A-Za-z]+)(?:,(\d+))?$")
_VARIANT_SUBSTITUTION = re.compile(r"^([A-Za-z]+),(\d+),([A-Za-z*?]+)$")


def _residue_sort_key(residue_name):
    """Sort key putting the nameless residue first, since `None < str` raises."""
    return (residue_name is not None, residue_name or "")


def _residue_name(amino_acid, residue):
    """`"S396"`, `"S"`, or `None` when the BEL element names no residue.

    `pmod(Ph)` says a protein is phosphorylated without saying where, and that is
    exactly a CellDesigner `<modificationResidue>` with no `name`: the writer
    omits the attribute when `residue.name is None`, and the reader reads it back
    as `None`, so the nameless residue survives a round trip. What must *not*
    happen is a modification with no **residue** at all -- the writer would emit
    `residue=""` and the reader looks that up unguarded -- so a residue object is
    always produced, named or not.
    """
    if amino_acid and residue:
        return f"{amino_acid}{residue}"
    if amino_acid:
        return amino_acid
    return None


def _read_protein_modification(protein_modification):
    """`(kind, payload)` for one `pmod`, before it is routed."""
    if protein_modification.namespace:
        # An ontology pmod -- `pmod(GO:"protein oxidation")` -- names a process,
        # not a residue state, so it cannot become a badge.
        return "structural_state", protein_modification.identifier
    state = PMOD_TYPE_TO_MODIFICATION_STATE.get(protein_modification.identifier)
    if state is None:
        raise RuntimeError(
            f"pmod type {protein_modification.identifier!r} has no CellDesigner "
            "modification state. Add it to PMOD_TYPE_TO_MODIFICATION_STATE."
        )
    return "residue_state", (
        _residue_name(protein_modification.amino_acid, protein_modification.residue),
        state,
    )


def _read_variant(variant):
    """`(kind, payload)` for one `var`, before it is routed.

    `?` and `p.?` occur verbatim in the data -- BEL saying "some variant,
    unspecified" -- and are kept: they are the only thing distinguishing that
    proteoform from the plain protein.
    """
    descriptor = re.sub(r"\s+", "", variant.descriptor or "")
    match = _VARIANT_MODIFICATION.match(descriptor)
    if match is not None:
        state = VARIANT_CODE_TO_MODIFICATION_STATE.get(match.group(1))
        if state is not None:
            return "residue_state", (
                _residue_name(match.group(2), match.group(3)),
                state,
            )
    match = _VARIANT_SUBSTITUTION.match(descriptor)
    if match is not None:
        return "structural_state", (
            f"{match.group(1)}{match.group(2)}{match.group(3)}"
        )
    return "structural_state", descriptor or None


def _read_fragment(fragment):
    """The structural state a `frag` contributes, or `None` when it says nothing.

    Nothing is ever invented: `frag("?")` -- BEL's "range unknown" -- says only
    "a fragment of", which the `TruncatedProtein` glyph already says, so it must
    not become a `StructuralState("?")` on top.
    """
    parts = []
    if fragment.start_stop and fragment.start_stop != "?":
        parts.append(fragment.start_stop)
    if fragment.descriptor:
        parts.append(fragment.descriptor)
    return ",".join(parts) if parts else None


def _name_suffix_part(kind, payload):
    """The text a modifier contributes to a species *name* when it cannot be drawn."""
    if kind == "residue_state":
        residue_name, state = payload
        return f"{state.value}{residue_name if residue_name is not None else ''}"
    return payload


def make_species_fields(elements):
    """`{element: fields}` -- what each element is, in CellDesigner vocabulary.

    `fields` is a plain dict, and it is what a species is built from *before any
    template exists*::

        species_class      the momapy class to instantiate
        template_class     its template's class, or None
        template_name      the template's name -- the bare identifier, so every
                           proteoform of one protein shares one template
        species_name       the glyph's label, identifier plus any modifier that
                           could not be drawn
        active             whether BEL asserts an activity of it
        compartment_name   its `loc()`, or None
        residue_states     ((residue name, ModificationState), ...) -- residue
                           *names*, not objects: the objects live on a template
                           that cannot be built until every proteoform has been
                           read
        structural_states  (value, ...)

    Covers `elements` **and** every element reached through their members, since
    a complex's subunits are species too. Pure and DB-free.

    An element that is the subject of some activity is `active`, exactly as the
    `Activity` of it is, so the two collapse into one species -- which is what
    keeps every upstream -> downstream path running through an activity and out
    of its subject (6889 of them at two hops in the AD KG, which three hops
    cannot recover). The price is that a protein carries an active border
    whenever BEL asserts any activity of it. The set is derived here, from the
    activities among the elements themselves.
    """
    all_elements = _collect_elements(elements)
    active_abundances = frozenset(
        element.abundance
        for element in all_elements
        if isinstance(element, momapy_bel.core.Activity)
    )
    species_fields = {}
    for element in all_elements:
        _add_species_fields(element, species_fields, active_abundances)
    return species_fields


def _member_elements(element):
    """The elements a species' subunits are built from, in a stable order.

    An activity stands for its subject, so it has its subject's members: the
    subunits of `act(complex(A,B))` are A and B. Its subject gets fields of its
    own through the delegation in :func:`_add_species_fields`.
    """
    if isinstance(element, momapy_bel.core.Activity):
        return _member_elements(element.abundance)
    return tuple(sorted(getattr(element, "members", ()), key=_member_sort_key))


def _member_sort_key(element):
    """A stable order for a complex's members, which are a frozenset."""
    return (
        type(element).__name__,
        element.namespace or "",
        element.identifier or "",
    )


def _collect_elements(elements):
    """`elements` plus every element reachable through their members, depth first."""
    collected = {}
    stack = list(elements)
    while stack:
        element = stack.pop()
        if element in collected:
            continue
        collected[element] = element
        stack.extend(_member_elements(element))
    return list(collected)


def _add_species_fields(element, species_fields, active_abundances):
    fields = species_fields.get(element)
    if fields is not None:
        return fields
    if isinstance(element, momapy_bel.core.Activity):
        # `act(p(X))` *is* X in an active state: it takes its subject's fields
        # whole, and the `ma()` code is dropped. Two activities of one subject
        # therefore give the same fields and collapse into one species --
        # accepted, and measured at 75 of the AD KG's 549.
        fields = _add_species_fields(
            element.abundance, species_fields, active_abundances
        )
    else:
        fields = _read_species_fields(element, species_fields, active_abundances)
    species_fields[element] = fields
    return fields


def _read_species_fields(element, species_fields, active_abundances):
    species_class, template_class = _resolve_cd_classes(element)
    modifiers = [
        _read_protein_modification(protein_modification)
        for protein_modification in getattr(element, "modifications", ())
    ] + [_read_variant(variant) for variant in getattr(element, "variants", ())]
    residue_states = []
    structural_states = []
    name_suffix_parts = []

    fragment = getattr(element, "fragment", None)
    if fragment is not None:
        # The class change happens even for a `frag("?")`: *that* is what says the
        # species is a fragment, and it says it whether or not the range is known.
        species_class, template_class = _TRUNCATED_CLASSES
        fragment_state = _read_fragment(fragment)
        if fragment_state is not None:
            structural_states.append(fragment_state)
    for kind, payload in modifiers:
        if payload is None:
            continue
        if kind == "residue_state" and issubclass(
            species_class, PROTEIN_SPECIES_CLASSES
        ):
            residue_states.append(payload)
        elif kind == "structural_state" and issubclass(
            species_class, PROTEIN_SPECIES_CLASSES + COMPLEX_SPECIES_CLASSES
        ):
            structural_states.append(payload)
        else:
            # A gene, an RNA, an abundance or a phenotype: no badge container, so
            # the modifier is folded into the species *name* -- never into the
            # template name, or one protein would declare N `<protein>` entries.
            name_suffix_parts.append(_name_suffix_part(kind, payload))

    if isinstance(element, momapy_bel.core.CompositeAbundance):
        # CellDesigner has no composite; the state says what the complex stands
        # for so a reader is not misled into thinking the members are bound.
        structural_states.append("composite")

    identifier = element.identifier
    if not identifier:
        identifier = ":".join(
            sorted(
                _add_species_fields(member, species_fields, active_abundances)[
                    "species_name"
                ]
                for member in _member_elements(element)
            )
        )
    species_name = identifier
    if name_suffix_parts:
        species_name = f"{identifier} [{'|'.join(sorted(set(name_suffix_parts)))}]"
    location = getattr(element, "location", None)
    return {
        "species_class": species_class,
        "template_class": template_class,
        "template_name": identifier if template_class is not None else None,
        "species_name": species_name,
        "active": element in active_abundances,
        "compartment_name": None if location is None else location.identifier,
        # Sorted through `_residue_sort_key`: a residue name may be `None` (a
        # `pmod(Ph)` naming no site), and `None` does not compare with a string.
        "residue_states": tuple(
            sorted(
                set(residue_states),
                key=lambda residue_state: (
                    _residue_sort_key(residue_state[0]),
                    residue_state[1].value,
                ),
            )
        ),
        "structural_states": tuple(sorted(set(structural_states))),
    }


def _resolve_cd_classes(element):
    if type(element) is momapy_bel.core.Abundance:
        if element.namespace in ABUNDANCE_MOLECULE_NAMESPACES:
            return momapy.celldesigner.SimpleMolecule, None
        return momapy.celldesigner.Unknown, None
    classes = BEL_CLASS_TO_CD_CLASSES.get(type(element))
    if classes is None:
        raise RuntimeError(
            f"BEL element {element!r} has class {type(element).__name__}, which has "
            "no CellDesigner counterpart. Add it to BEL_CLASS_TO_CD_CLASSES."
        )
    return classes


# ---------------------------------------------------------------------------
# templates and compartments: the shared objects, interned by value
# ---------------------------------------------------------------------------


def make_templates(species_fields):
    """`{(template class, template name): template}` -- one object per key.

    That is invariant (a) within the export, and invariant (b) rests on it: a
    template's residues are the union over **every** proteoform of that protein,
    which is why they cannot be built one species at a time.

    Identity with the *stored* collections' templates is not this function's job:
    the export is standalone, and hash integration with a seeded
    `object_key_to_node` makes two value-equal elements one database node across
    save calls.
    """
    residue_names = {}
    for fields in species_fields.values():
        template_class = fields["template_class"]
        if template_class is None:
            continue
        names = residue_names.setdefault(
            (template_class, fields["template_name"]), set()
        )
        if template_class in TEMPLATE_RESIDUE_FIELDS:
            names.update(residue_name for residue_name, _ in fields["residue_states"])

    templates = {}
    for key in sorted(residue_names, key=lambda key: (key[0].__name__, key[1])):
        template_class, template_name = key
        field_name, residue_class = TEMPLATE_RESIDUE_FIELDS[template_class]
        # One residue object per (name, order) in one place: two value-equal
        # residue objects would collapse in the template's frozenset *after*
        # `renumber_ids` gave them different ids, where nothing is left to
        # notice -- and the modification would then name a residue the
        # `<listOfModificationResidues>` does not declare.
        template_residues = frozenset(
            residue_class(name=residue_name, order=order)
            for order, residue_name in enumerate(
                sorted(residue_names[key], key=_residue_sort_key)
            )
        )
        templates[key] = template_class(
            name=template_name, **{field_name: template_residues}
        )
    return templates


def make_compartments(species_fields, root_compartment_name):
    """`({location name: Compartment}, root)` -- one object per location.

    Invariant (d). The root is undrawn, like the stored maps' `default`: it
    exists to give every BEL compartment an `outside`, which is what makes
    `pd2af.utils._build_dot_graph` attach their clusters at all, and what keeps a
    BEL `nucleus` value-distinct from a stored one.
    """
    compartment_root = momapy.celldesigner.Compartment(name=root_compartment_name)
    compartments = {
        name: momapy.celldesigner.Compartment(name=name, outside=compartment_root)
        for name in sorted(
            {
                fields["compartment_name"]
                for fields in species_fields.values()
                if fields["compartment_name"] is not None
            }
        )
    }
    return compartments, compartment_root


# ---------------------------------------------------------------------------
# the species
# ---------------------------------------------------------------------------


def make_species(
    element,
    species_fields,
    templates,
    compartments,
    is_subunit=False,
    record=None,
):
    """One frozen momapy species for one BEL element. Always a **fresh object tree**.

    Invariant (c): the writer indexes complex subunits by `id()` and skips any
    top-level species whose `id()` is in that index, so memoising a species by
    element -- the obvious way to get value-collapse -- would make the own glyph
    of every protein that is also a complex member disappear, and reroute every
    modulation targeting it to the enclosing complex. Interning happens only at
    the top level, on the finished object, in `commute_dm.bel_export`.

    A subunit keeps no compartment: CellDesigner puts an included species inside
    its complex, not inside a compartment box, and only top-level species
    contribute to `model.compartments`.

    When `record` is a list, `(element, species)` is appended for **every**
    species built, subunits included. That is what lets a subunit be annotated: a
    subunit is a per-occurrence object, and 16 of the interface's UniProt
    identifiers are carried only by complex members.
    """
    fields = species_fields[element]
    species_class = fields["species_class"]
    template_class = fields["template_class"]
    arguments = {"name": fields["species_name"], "active": fields["active"]}
    if not is_subunit and fields["compartment_name"] is not None:
        arguments["compartment"] = compartments[fields["compartment_name"]]
    if template_class is not None:
        template = templates[(template_class, fields["template_name"])]
        arguments["template"] = template
        arguments["modifications"] = frozenset()
        if template_class in TEMPLATE_RESIDUE_FIELDS:
            field_name, _ = TEMPLATE_RESIDUE_FIELDS[template_class]
            # Invariant (b): the residue objects come from the template the
            # species itself carries, never from a value-equal copy.
            residue_by_name = {
                residue.name: residue for residue in getattr(template, field_name)
            }
            arguments["modifications"] = frozenset(
                momapy.celldesigner.Modification(
                    residue=residue_by_name[residue_name], state=state
                )
                for residue_name, state in fields["residue_states"]
            )
    if issubclass(species_class, PROTEIN_SPECIES_CLASSES + COMPLEX_SPECIES_CLASSES):
        arguments["structural_states"] = frozenset(
            momapy.celldesigner.StructuralState(value=value)
            for value in fields["structural_states"]
        )
    if issubclass(species_class, COMPLEX_SPECIES_CLASSES):
        arguments["subunits"] = frozenset(
            make_species(
                member,
                species_fields,
                templates,
                compartments,
                is_subunit=True,
                record=record,
            )
            for member in _member_elements(element)
        )
    species = species_class(**arguments)
    if record is not None:
        record.append((element, species))
    return species


# ---------------------------------------------------------------------------
# the glyphs
# ---------------------------------------------------------------------------

# Label measuring and wrapping live in `commute_dm.submaps`: `commute_dm.core`
# needs `make_fitted_synthetic_layout` for every sub-map's synthetic central
# node, in both pairings, and this module is export-only. They are imported back
# under their old names; `submaps` imports nothing of this module, so the
# dependency is acyclic.
_LABEL_MAX_WIDTH = commute_dm.submaps._LABEL_MAX_WIDTH
_LABEL_PADDING = commute_dm.submaps._LABEL_PADDING
_MEASURING_POSITION = commute_dm.submaps._MEASURING_POSITION
_measure = commute_dm.submaps._measure
_wrap_label_text = commute_dm.submaps._wrap_label_text
make_fitted_synthetic_layout = commute_dm.submaps.make_fitted_synthetic_layout

# Subunit stacking inside a complex, and the band its own label sits in.
_SUBUNIT_XSEP = 12.0
_SUBUNIT_YSEP = 12.0
_COMPLEX_LABEL_BAND = 22.0
# A badge straddles its species' border, so a decorated glyph is grown to keep
# it off its neighbours -- `make_auto_layout` spaces nodes by the size the
# layout element carries.
_BADGE_MARGIN = 16.0
# Horizontal gap between two glyphs' throwaway synthetic positions. Only its
# being positive matters -- see :func:`make_species_layouts`.
_GLYPH_GAP = 20.0

# The `*ActiveLayout` sibling the reader appends to an active species' glyph.
# Emitting it forward makes a second write of a read-back map stable.
_LAYOUT_CLASS_TO_ACTIVE_LAYOUT_CLASS = {
    momapy.celldesigner.GenericProteinLayout: (
        momapy.celldesigner.GenericProteinActiveLayout
    ),
    momapy.celldesigner.TruncatedProteinLayout: (
        momapy.celldesigner.TruncatedProteinActiveLayout
    ),
    momapy.celldesigner.ComplexLayout: momapy.celldesigner.ComplexActiveLayout,
    momapy.celldesigner.GeneLayout: momapy.celldesigner.GeneActiveLayout,
    momapy.celldesigner.RNALayout: momapy.celldesigner.RNAActiveLayout,
    momapy.celldesigner.SimpleMoleculeLayout: (
        momapy.celldesigner.SimpleMoleculeActiveLayout
    ),
    momapy.celldesigner.UnknownLayout: momapy.celldesigner.UnknownActiveLayout,
    momapy.celldesigner.PhenotypeLayout: momapy.celldesigner.PhenotypeActiveLayout,
}


@dataclasses.dataclass(frozen=True)
class _Box:
    """The measured size of one species' glyph and, recursively, its subunits."""

    species: object
    width: float
    height: float
    label_text: str
    children: tuple = ()


def _species_sort_key(species):
    """A stable order for a complex's subunits, which are a frozenset."""
    return (
        species.name or "",
        type(species).__name__,
        tuple(sorted(s.value or "" for s in getattr(species, "structural_states", ()))),
        tuple(
            sorted(
                (m.residue.name or "", m.state.value if m.state else "")
                for m in getattr(species, "modifications", ())
            )
        ),
    )


def _has_decorations(species):
    return bool(
        getattr(species, "modifications", ())
        or getattr(species, "structural_states", ())
    )


def _measure_species(species):
    """The :class:`_Box` of `species`, measuring subunits before their complex.

    A complex is sized to its stacked subunits, a leaf to its wrapped label, and
    both to at least their layout class's default. A decorated glyph is grown by
    `_BADGE_MARGIN` on both axes so its badges, which straddle the border, do
    not land on a neighbour.
    """
    default = pd2af.celldesigner.building_layout.make_synthetic_layout(species, 0)
    label_text, label_width, label_height = _wrap_label_text(
        default.label.text if default.label is not None else ""
    )
    children = tuple(
        _measure_species(subunit)
        for subunit in sorted(
            getattr(species, "subunits", ()) or (), key=_species_sort_key
        )
    )
    if children:
        inner_width = max(child.width for child in children)
        inner_height = sum(child.height for child in children) + _SUBUNIT_YSEP * (
            len(children) - 1
        )
        width = max(
            default.width, inner_width + 2 * _SUBUNIT_XSEP, label_width + _LABEL_PADDING
        )
        height = max(
            default.height,
            inner_height + 2 * _SUBUNIT_YSEP + _COMPLEX_LABEL_BAND,
        )
    else:
        width = max(default.width, label_width + _LABEL_PADDING)
        height = max(default.height, label_height + _LABEL_PADDING)
    if _has_decorations(species):
        width += _BADGE_MARGIN
        height += _BADGE_MARGIN
    return _Box(
        species=species,
        width=width,
        height=height,
        label_text=label_text,
        children=children,
    )


def _make_badge_label(badge, text):
    return dataclasses.replace(
        badge,
        label=momapy.core.layout.TextLayout(
            text=text,
            font_size=momapy.celldesigner.DEFAULT_MODIFICATION_FONT_SIZE,
            font_family=momapy.drawing.DEFAULT_FONT_FAMILY,
            fill=momapy.coloring.black,
            stroke=momapy.drawing.NoneValue,
            position=badge.label_center(),
        ),
    )


def _badge_position(node, angle):
    """The reader's exact forward transform from a CellDesigner angle to a point.

    Kept identical so that `_writing.compute_cd_angle` inverts it and the angle
    written out is the one intended.
    """
    point = momapy.geometry.Point(
        node.width * math.cos(angle), node.height * math.sin(angle)
    )
    position = node.own_angle(math.atan2(point.y, point.x), unit="radians")
    return position if position is not None else node.center()


def _make_decoration_layouts(species, glyph, mapping):
    """The `ModificationLayout` / `StructuralStateLayout` children of `glyph`.

    **Badge angles are residue-derived, not species-local.** CellDesigner stores
    the angle on the *template* and `_writing.find_residue_angle` takes it from
    the first species using that template with a modification on that residue,
    so a species-local index would move every other proteoform's badge on
    read-back. The angle is `2*pi * residue.order / <number of the template's
    residues>`, which every proteoform of one protein agrees on.

    Structural states carry no angle in the XML, so they go on a disjoint arc
    with a species-local index.
    """
    layout_elements = []
    modifications = sorted(
        getattr(species, "modifications", ()) or (),
        key=lambda modification: (
            modification.residue.order,
            modification.state.value if modification.state else "",
        ),
    )
    if modifications:
        field_name, _ = TEMPLATE_RESIDUE_FIELDS[type(species.template)]
        n_residues = len(getattr(species.template, field_name)) or 1
        for modification in modifications:
            angle = 2 * math.pi * modification.residue.order / n_residues
            badge = momapy.celldesigner.ModificationLayout(
                position=_badge_position(glyph, angle)
            )
            badge = _make_badge_label(
                badge, modification.state.value if modification.state else ""
            )
            mapping.add_mapping(badge, modification)
            layout_elements.append(badge)
    structural_states = sorted(
        getattr(species, "structural_states", ()) or (),
        key=lambda structural_state: structural_state.value or "",
    )
    for index, structural_state in enumerate(structural_states):
        angle = 2 * math.pi * (index + 0.5) / len(structural_states)
        badge = momapy.celldesigner.StructuralStateLayout(
            position=_badge_position(glyph, angle)
        )
        badge = _make_badge_label(badge, structural_state.value or "")
        mapping.add_mapping(badge, structural_state)
        layout_elements.append(badge)
    return layout_elements


def _make_active_layouts(species, glyph):
    if not species.active:
        return []
    active_class = _LAYOUT_CLASS_TO_ACTIVE_LAYOUT_CLASS.get(type(glyph))
    if active_class is None:
        raise RuntimeError(
            f"{type(glyph).__name__} has no active-border layout class; add it to "
            "_LAYOUT_CLASS_TO_ACTIVE_LAYOUT_CLASS."
        )
    return [
        active_class(
            position=glyph.position,
            width=glyph.width + 2 * momapy.celldesigner.DEFAULT_ACTIVE_XSEP,
            height=glyph.height + 2 * momapy.celldesigner.DEFAULT_ACTIVE_YSEP,
        )
    ]


def make_species_layouts(species_list, mapping):
    """One glyph per species, with subunits nested and badges attached.

    Two passes per glyph -- measure, then build -- so that no layout element is
    `dataclasses.replace`d after it has entered `mapping`: the mapping's forward
    dict is equality-keyed, and a replaced key would be a *different* key.

    Registering `mapping.add_mapping(subunit_layout, subunit)` for every subunit
    is the load-bearing part: the writer emits an included species' alias only
    when `get_child_layout_elements(subunit, complex)` returns a layout that is
    **also** in `complex_layout.layout_elements`, so both the nesting and the
    mapping entry are required. (This is what refutes the old "a subunit's
    sub-glyph does not survive a CellDesigner round trip" comment.)

    **Why the glyphs are laid out here rather than one at a time.** The positions
    are throwaway -- `make_auto_layout` recomputes every one of them -- but they
    are not arbitrary: a badge is a 16x16 white square whose only distinguishing
    content is its position and its one-character label, so two badges of two
    species at nearby positions are *value-equal*, and the second
    `add_mapping(badge, modification)` then evicts the first from the
    equality-keyed forward dict. The species are therefore spread along x with a
    gap, so that the x-intervals `[center - width/2, center + width/2]` of any
    two glyphs are disjoint and no two decorations of different species can
    coincide. (Measured: with pd2af's synthetic 1pt spacing, one AD modification
    silently lost its mapping entry.)

    Nested geometry *does* survive `make_auto_layout`: `_build_dot_graph` makes
    dot nodes only from top-level layout elements and takes their size from the
    layout element, and `_reposition_from_dot` translates the whole subtree.
    """
    boxes = [_measure_species(species) for species in species_list]
    layout_elements = []
    cursor = 0.0
    for box in boxes:
        cursor += box.width / 2
        layout_elements.append(
            _build_species_layout(box, momapy.geometry.Point(cursor, 0.0), mapping)
        )
        cursor += box.width / 2 + _GLYPH_GAP
    return layout_elements


def _build_species_layout(box, center, mapping):
    species = box.species
    children = []
    if box.children:
        content_height = sum(child.height for child in box.children) + _SUBUNIT_YSEP * (
            len(box.children) - 1
        )
        # The complex's own label sits in a band at the bottom, so the subunit
        # stack is centred above it.
        top = center.y - _COMPLEX_LABEL_BAND / 2 - content_height / 2
        for child_box in box.children:
            child_center = momapy.geometry.Point(
                center.x, top + child_box.height / 2
            )
            child_layout = _build_species_layout(
                child_box, child_center, mapping
            )
            mapping.add_mapping(child_layout, child_box.species)
            children.append(child_layout)
            top += child_box.height + _SUBUNIT_YSEP
    glyph = pd2af.celldesigner.building_layout.make_synthetic_layout(species, 0)
    glyph = dataclasses.replace(
        glyph, position=center, width=box.width, height=box.height
    )
    if glyph.label is not None:
        glyph = dataclasses.replace(
            glyph,
            label=dataclasses.replace(
                glyph.label, text=box.label_text, position=glyph.label_center()
            ),
        )
    # Badges and the active border are measured against the *childless* glyph,
    # which is right: `own_angle` and `center` ignore `layout_elements`.
    decorations = _make_decoration_layouts(species, glyph, mapping)
    actives = _make_active_layouts(species, glyph)
    return dataclasses.replace(
        glyph, layout_elements=tuple(children) + tuple(decorations) + tuple(actives)
    )


def make_compartment_layout(compartment):
    """A drawn box for one `loc()` compartment.

    Geometry is throwaway: `make_auto_layout`'s
    `_apply_dot_cluster_bounding_boxes_to_compartments` copies each dot cluster's
    bounding box onto its compartment layout. The **BEL root compartment gets no
    layout**, exactly like the stored maps' `default` root -- it is a pure
    grouping parent, so no box is drawn for it and the visual result is
    top-level `loc()` boxes.
    """
    label = momapy.core.layout.TextLayout(
        text=compartment.name or "",
        position=_MEASURING_POSITION,
    )
    return momapy.celldesigner.RectangleCompartmentLayout(
        position=_MEASURING_POSITION,
        width=1.0,
        height=1.0,
        label=label,
    )
