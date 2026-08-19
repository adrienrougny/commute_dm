"""BEL terms as real activity-flow CellDesigner content.

`commute_dm.bel_export` turns a BEL knowledge graph's influence-graph
projection into an in-memory `CellDesignerMap`. This module is the half that
knows what a BEL *term* is: it reads the term graph out of Neo4j, describes each
term in CellDesigner vocabulary, and builds the frozen momapy species and glyphs.
`bel_export` imports `bel_terms`, never the reverse.

**The structure is fully relational.** Every sub-term of a BEL term -- a `pmod`,
a `var`, a `frag`, a `loc`, each complex member, an activity's subject -- is its
own Neo4j node reached by a typed `HAS__*` edge, and the entity nodes carry
`name` / `namespace` properties. Nothing here parses a BEL string.

What a term becomes::

    p(HGNC:"GDE1")                        GenericProtein("GDE1")
    p(HGNC:"MAPT",pmod(Ph,Y,18))          + Modification(residue Y18, PHOSPHORYLATED)
    p(HGNC:"MAPT",var("P, T, 212"))       + Modification(residue T212, PHOSPHORYLATED)
    p(MGI:"App",var("K,670,N"))           + StructuralState("K670N")
    p(HGNC:"APP",frag("672_713"))         TruncatedProtein("APP") + StructuralState
    p(HGNC:"CLU",loc(GO:"intracellular")) compartment=Compartment("intracellular")
    act(p(HGNC:"CSNK2A1"),ma(kin))        GenericProtein("CSNK2A1", active=True)
    complex(p(A),p(B))                    Complex("A:B", subunits={A, B})
    composite(p(A),p(B))                  Complex + StructuralState("composite")
    a(CHEBI:"dopamine")                   SimpleMolecule("dopamine")
    a(CONSO:"Tau aggregates")             Unknown("Tau aggregates")
    bp(GO:...) / path(MESH:...)           Phenotype

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
:func:`make_build_context`. Identity with the *stored* collections' templates is
not established here: the export is standalone, and the save merges value-equal
elements into one database node (hash integration with a seeded
`object_key_to_node`, see `2_10`).

(b) **`Modification.residue` *is* an object inside its species' template's
residue container.** The writer emits the modification's residue id and the
template's residue ids independently; after `renumber_ids` two value-equal
residue objects have two different ids and the modification names a residue that
is not in the list, which the reader looks up unguarded. Hence a template's
residues are the union over *every* proteoform of that protein, and the residue
lookup is re-derived from the *interned* template.

(c) **A subunit object is never the same object as a top-level species.** The
writer keys `build_subunit_to_complex` by `id(subunit)` and skips any species
whose `id()` is in that index. `p(HGNC:"APP")` is both a projected node and a
member of several complexes, so memoising species by node id would make APP's own
glyph disappear. :func:`make_species` therefore always builds a **fresh object
tree**; interning happens only at the top level, in `bel_export`.

This one is an **in-memory property of the export only, and does not survive the
save**: a top-level `p(X)` and the same `p(X)` inside a complex are value-equal,
so hash integration stores them as one node and hydration returns one object.
A source map built from the stored collection therefore has hundreds of subunits
that *are* their top-level species (651 in AD), and that is expected --
`submaps.make_submap_from_model_elements` promotes such a species to its own
object per sub-map, before the invariant checks run. Never run
`check_identity_invariants` on a source map; it is a pre-write check for
sub-maps.

(d) **One compartment object per `(collection, location)`, hanging off an undrawn
per-collection root.** Same collapse as (a), one level down -- and the root
parent is not decorative: `pd2af.utils._build_dot_graph` only attaches a
compartment's dot cluster to the graph when `compartment.outside is not None`
(there is no `else`), so a drawn compartment with no `outside` gets an orphan
cluster, its species never reach graphviz, and they keep their throwaway
positions. The root also keeps BEL compartments structurally distinct from the
stored CellDesigner ones, 18 of whose names a BEL `loc()` shares: a BEL
`nucleus` has `outside=<BEL root>`, a stored one `outside=default`, so they can
never be value-equal and neither can the species inside them.
"""

import dataclasses
import math
import re

import momapy.celldesigner
import momapy.coloring
import momapy.core.layout
import momapy.drawing
import momapy.geometry
import pd2af.celldesigner.building_layout

import commute_dm.submaps


# ---------------------------------------------------------------------------
# the term graph
# ---------------------------------------------------------------------------

# The concrete BEL node classes, most specific first. pylpg labels a node with
# its own class *and* every ancestor (`['GeneticFlow', 'Protein',
# 'BELModelElement']`), so the concrete class is resolved here rather than by
# taking `labels(node)[0]`, which is order-dependent. Exactly one of these must
# match, so an unexpected label combination fails loudly.
BEL_TERM_NODE_CLASSES = (
    "MicroRna",
    "Rna",
    "Gene",
    "Protein",
    "Complex",
    "Composite",
    "Abundance",
    "Activity",
    "BiologicalProcess",
    "Pathology",
    "Population",
    "Variant",
    "ProteinModification",
    "Fragment",
    "Gmod",
    "ToLocation",
    "FromLocation",
    "Location",
    "Translocation",
    "Degradation",
    "CellSecretion",
    "CellSurfaceExpression",
    "Reaction",
    "Reactants",
    "Products",
    "List",
)

# Nodes a `BELModel` holds that are not terms at all. `Subgraph` is the KG's own
# grouping plumbing (`HAS_NODE` reaches it like any other node).
IGNORED_BEL_NODE_CLASSES = frozenset({"Subgraph"})

# A complex's or composite's members, and an activity's subject: the sub-terms
# that are entities in their own right.
MEMBER_EDGE_TYPES = frozenset(
    {
        "HAS__PROTEIN",
        "HAS__ABUNDANCE",
        "HAS__COMPLEX",
        "HAS__COMPOSITE",
        "HAS__GENE",
        "HAS__RNA",
        "HAS__MICRO_RNA",
    }
)
# The sub-terms that qualify an entity rather than being one.
MODIFIER_EDGE_TYPES = frozenset(
    {
        "HAS__PMOD",
        "HAS__VARIANT",
        "HAS__FRAGMENT",
        "HAS__LOCATION",
        "HAS__GMOD",
    }
)
# `HAS__*` edges of BEL *process* terms (translocations, reactions, degradations,
# lists). None of these terms is projected, so none is ever described -- but they
# are enumerated so that :func:`load_bel_terms` can raise on an edge type it has
# never seen instead of silently building a term with a missing part.
OTHER_STRUCTURAL_EDGE_TYPES = frozenset(
    {
        "HAS__TO_LOCATION",
        "HAS__FROM_LOCATION",
        "HAS__REACTANTS",
        "HAS__PRODUCTS",
    }
)
STRUCTURAL_EDGE_TYPES = (
    MEMBER_EDGE_TYPES | MODIFIER_EDGE_TYPES | OTHER_STRUCTURAL_EDGE_TYPES
)

# Every `:BELModel` node of the collections, with the properties the term
# description reads. `default` (the activity code), `type` (the pmod type),
# `position` and `range` are Cypher reserved-ish words and must be backticked.
# `DISTINCT` because two collection entries can reach the same model.
_TERM_NODES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(node)
WHERE collection.name IN $collection_names
RETURN DISTINCT elementId(node) AS node_id,
       collection.name AS collection_name,
       labels(node) AS node_labels,
       node.bel AS bel,
       node.name AS name,
       node.namespace AS namespace,
       node.`default` AS activity_code,
       node.hgvs AS hgvs,
       node.`type` AS modification_type,
       node.amino_acid AS amino_acid,
       node.`position` AS position,
       node.`range` AS range,
       node.descriptor AS descriptor
"""

# The structural (sub-term) edges. The **double** underscore is deliberate: it
# excludes the sparse single-underscore statement route (`HAS_MEMBERS`,
# `HAS_COMPONENTS`, `HAS_MODIFIED_GENE`, ...) and the `HAS_NODE` / `HAS_OBJ`
# plumbing.
_TERM_EDGES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(container)
MATCH (container)-[part]->(node)
WHERE collection.name IN $collection_names
    AND type(part) STARTS WITH 'HAS__'
RETURN DISTINCT elementId(container) AS container_node_id,
       elementId(node) AS node_id,
       type(part) AS edge_type
"""


@dataclasses.dataclass(eq=False)
class BelTerm:
    """One BEL term node and the sub-terms it points at.

    `eq=False` on purpose: the term graph is recursive and two content-equal
    terms are two different BEL nodes, so identity is what we want -- it is also
    what `describe_term`'s memo and `make_build_context` are keyed on.
    """

    node_id: str
    collection_name: str
    node_class: str
    bel: str | None = None
    name: str | None = None
    namespace: str | None = None
    activity_code: str | None = None
    hgvs: str | None = None
    modification_type: str | None = None
    amino_acid: str | None = None
    position: str | None = None
    range: str | None = None
    descriptor: str | None = None
    parts: list = dataclasses.field(default_factory=list)

    def members(self):
        """The sub-terms that are entities: complex members, activity subject."""
        return [part for edge_type, part in self.parts if edge_type in MEMBER_EDGE_TYPES]

    def modifiers(self):
        """The `(edge_type, sub-term)` pairs that qualify this term."""
        return [
            (edge_type, part)
            for edge_type, part in self.parts
            if edge_type in MODIFIER_EDGE_TYPES
        ]


def _resolve_node_class(node_id, node_labels):
    matched = [label for label in BEL_TERM_NODE_CLASSES if label in node_labels]
    if len(matched) == 1:
        return matched[0]
    if not matched:
        raise RuntimeError(
            f"BEL node {node_id} has labels {sorted(node_labels)}, none of which is a "
            "known term class. Add it to BEL_TERM_NODE_CLASSES (or to "
            "IGNORED_BEL_NODE_CLASSES if it is not a term)."
        )
    raise RuntimeError(
        f"BEL node {node_id} has labels {sorted(node_labels)}, matching several term "
        f"classes {matched}. BEL_TERM_NODE_CLASSES must be ordered most specific "
        "first and must not contain two classes one node can carry at once."
    )


def load_bel_terms(session, collection_names):
    """`{node_id: BelTerm}` for every term node of the given BEL collections.

    Two queries, ~5000 nodes and ~6000 edges, well under a second. Node ids are
    Neo4j `elementId`s -- the same keys `commute_dm.submaps.Influences` and
    momapy_kb's `node_id_to_object` cache use, which is what lets the projection
    and the term graph be joined without a second index.

    Raises on a node class or a `HAS__*` edge type this module does not know:
    silently skipping is how the old influence-graph code rotted into a smaller
    graph than it claimed.
    """
    parameters = {"collection_names": list(collection_names)}
    terms = {}
    for row in session.execute_query(_TERM_NODES_QUERY, parameters):
        node_labels = row["node_labels"]
        if any(label in IGNORED_BEL_NODE_CLASSES for label in node_labels):
            continue
        node_id = row["node_id"]
        terms[node_id] = BelTerm(
            node_id=node_id,
            collection_name=row["collection_name"],
            node_class=_resolve_node_class(node_id, node_labels),
            bel=row["bel"],
            name=row["name"],
            namespace=row["namespace"],
            activity_code=row["activity_code"],
            hgvs=row["hgvs"],
            modification_type=row["modification_type"],
            amino_acid=row["amino_acid"],
            position=None if row["position"] is None else str(row["position"]),
            range=row["range"],
            descriptor=row["descriptor"],
        )
    for row in session.execute_query(_TERM_EDGES_QUERY, parameters):
        edge_type = row["edge_type"]
        if edge_type not in STRUCTURAL_EDGE_TYPES:
            raise RuntimeError(
                f"BEL node {row['container_node_id']} has a {edge_type} sub-term edge, "
                "which is not classified. Add it to MEMBER_EDGE_TYPES, to "
                "MODIFIER_EDGE_TYPES or to OTHER_STRUCTURAL_EDGE_TYPES."
            )
        container = terms.get(row["container_node_id"])
        part = terms.get(row["node_id"])
        if container is None or part is None:  # an ignored (non-term) endpoint
            continue
        container.parts.append((edge_type, part))
    # Sorted so that a given knowledge graph always yields the same subunit
    # order, the same synthetic positions and the same `renumber_ids` numbering.
    for term in terms.values():
        term.parts.sort(key=lambda pair: (pair[0], pair[1].name or "", pair[1].node_id))
    return terms


# ---------------------------------------------------------------------------
# describing a term in CellDesigner vocabulary
# ---------------------------------------------------------------------------

# The BEL entity classes that have a fixed CellDesigner counterpart, mapped to
# `(species class, template class)`. `Activity` is absent -- it has no class of
# its own and takes its subject's -- and so is `Abundance`, which splits on its
# `namespace` property (see below).
#
# This table is **not** cosmetic any more: it selects the template class too, so
# an entry decides both the glyph and whether the species declares a
# `<protein>` / `<gene>` / `<rna>`.
#
# A microRNA is RNA; `AntisenseRNA` would assert a strand the KG does not give.
# CellDesigner has no "composite", and a set of abundances is closest to a
# complex; it has no disease class either, and a pathology is a state of the
# system, so both `bp()` and `path()` are phenotypes.
BEL_CLASS_TO_CD_CLASSES = {
    "Protein": (
        momapy.celldesigner.GenericProtein,
        momapy.celldesigner.GenericProteinTemplate,
    ),
    "Gene": (momapy.celldesigner.Gene, momapy.celldesigner.GeneTemplate),
    "Rna": (momapy.celldesigner.RNA, momapy.celldesigner.RNATemplate),
    "MicroRna": (momapy.celldesigner.RNA, momapy.celldesigner.RNATemplate),
    "Complex": (momapy.celldesigner.Complex, None),
    "Composite": (momapy.celldesigner.Complex, None),
    "BiologicalProcess": (momapy.celldesigner.Phenotype, None),
    "Pathology": (momapy.celldesigner.Phenotype, None),
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
GENE_RNA_SPECIES_CLASSES = (momapy.celldesigner.Gene, momapy.celldesigner.RNA)

# Where a template keeps the residues a `Modification` can name, and which class
# those residues are. A `GeneTemplate` / `RNATemplate` keeps them in `regions` as
# `ModificationSite`s instead, which the writer emits differently -- gene and RNA
# badges are deferred, so enabling them later is one entry here plus dropping
# `GENE_RNA_SPECIES_CLASSES` from the name-suffix branch.
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

# A BEL `pmod` type code to the CellDesigner state its badge shows. A type
# outside this table is an ontology pmod (`pmod(GO:"protein oxidation")`), which
# names a process rather than a residue state and becomes a structural state.
PMOD_TYPE_TO_MODIFICATION_STATE = {
    "pho": momapy.celldesigner.ModificationState.PHOSPHORYLATED,
    "ubi": momapy.celldesigner.ModificationState.UBIQUITINATED,
    "ace": momapy.celldesigner.ModificationState.ACETYLATED,
    "me0": momapy.celldesigner.ModificationState.METHYLATED,
    "gly": momapy.celldesigner.ModificationState.GLYCOSYLATED,
    "ogl": momapy.celldesigner.ModificationState.GLYCOSYLATED,
}

# BEL `var()` strings are free text in practice. Two shapes are recognised, in
# this order, after whitespace is stripped (so `"P, S, 9"` and `"P,S,9"` are the
# same variant): a phospho-form, which is really a `pmod` written as a variant,
# and an amino-acid substitution. Anything else is kept verbatim as a structural
# state -- `misfolded`, `p.D614G`, `c.863G>A`, `del`.
_VARIANT_PHOSPHO = re.compile(r"^(?:P|Ph),([A-Za-z]+)(?:,(\d+))?$")
_VARIANT_SUBSTITUTION = re.compile(r"^([A-Za-z]+),(\d+),([A-Za-z*?]+)$")


@dataclasses.dataclass(frozen=True)
class TermDescription:
    """What one BEL term is, in CellDesigner vocabulary. Pure data, no objects.

    `members` are the `BelTerm`s of a complex's or composite's subunits, kept as
    terms rather than as descriptions so :func:`make_species` can build a fresh
    object tree per occurrence -- invariant (c).
    """

    species_class: type
    template_class: type | None
    template_name: str | None
    species_name: str
    active: bool = False
    compartment_name: str | None = None
    modifications: tuple = ()  # ((residue_name, ModificationState | None), ...)
    structural_states: tuple = ()
    members: tuple = ()


def _residue_sort_key(residue_name):
    """Sort key putting the nameless residue first, since `None < str` raises."""
    return (residue_name is not None, residue_name or "")


def _residue_name(amino_acid, position):
    """`"S396"`, `"S"`, or `None` when the BEL term names no residue.

    `pmod(Ph)` says a protein is phosphorylated without saying where, and that is
    exactly a CellDesigner `<modificationResidue>` with no `name`: the writer
    omits the attribute when `residue.name is None`, and the reader reads it back
    as `None`, so the nameless residue survives a round trip. What must *not*
    happen is a modification with no **residue** at all -- the writer would emit
    `residue=""` and the reader looks that up unguarded -- so a residue object is
    always produced, named or not.
    """
    if amino_acid and position:
        return f"{amino_acid}{position}"
    if amino_acid:
        return amino_acid
    return None


def _describe_modifier(edge_type, part):
    """`(kind, payload)` for one modifier sub-term, before it is routed.

    `kind` is one of `"modification"` (payload `(residue name, state)`),
    `"structural_state"` (payload a string) or `"fragment"` (payload a string).
    Routing to the species' actual capabilities happens in :func:`describe_term`.

    A **`None` payload means the sub-term carries no information to draw**, and
    the caller records nothing for it. Nothing is ever invented: a `frag("?")`
    says only "a fragment of", which the `TruncatedProtein` glyph already says,
    so it must not become a `StructuralState("?")`. A `var("?")` is different --
    there `"?"` *is* the value BEL states, and it is the only thing telling that
    proteoform apart from the plain protein, so it is kept.
    """
    if edge_type == "HAS__PMOD":
        state = PMOD_TYPE_TO_MODIFICATION_STATE.get(part.modification_type)
        if state is not None:
            return "modification", (
                _residue_name(part.amino_acid, part.position),
                state,
            )
        if part.modification_type is not None:
            raise RuntimeError(
                f"BEL node {part.node_id} ({part.bel!r}) has pmod type "
                f"{part.modification_type!r}, which has no CellDesigner modification "
                "state. Add it to PMOD_TYPE_TO_MODIFICATION_STATE."
            )
        # An ontology pmod -- `pmod(GO:"protein oxidation")` -- names a process,
        # not a residue state, so it cannot become a badge.
        return "structural_state", part.name or part.bel or None
    if edge_type == "HAS__VARIANT":
        variant = re.sub(r"\s+", "", part.hgvs or "")
        match = _VARIANT_PHOSPHO.match(variant)
        if match is not None:
            return "modification", (
                _residue_name(match.group(1), match.group(2)),
                momapy.celldesigner.ModificationState.PHOSPHORYLATED,
            )
        match = _VARIANT_SUBSTITUTION.match(variant)
        if match is not None:
            return "structural_state", (
                f"{match.group(1)}{match.group(2)}{match.group(3)}"
            )
        # `?` and `p.?` occur verbatim in the data -- BEL saying "some variant,
        # unspecified" -- and are kept: they are the only thing distinguishing
        # that proteoform from the plain protein.
        return "structural_state", variant or None
    if edge_type == "HAS__FRAGMENT":
        # `?` is BEL's "range unknown", not a range: it says nothing the
        # `TruncatedProtein` glyph does not already say.
        parts = []
        if part.range and part.range != "?":
            parts.append(part.range)
        if part.descriptor:
            parts.append(part.descriptor)
        return "fragment", ",".join(parts) if parts else None
    if edge_type == "HAS__GMOD":
        return "structural_state", part.name or part.bel or None
    raise RuntimeError(
        f"BEL node {part.node_id} ({part.bel!r}) is reached by a {edge_type} modifier "
        "edge that _describe_modifier does not handle."
    )


def _name_suffix_part(kind, payload):
    """The text a modifier contributes to a species *name* when it cannot be drawn."""
    if kind == "modification":
        residue_name, state = payload
        return f"{state.value}{residue_name if residue_name is not None else ''}"
    return payload


def collect_activity_subject_term_ids(terms):
    """`{id(term)}` for every term that is the subject of some activity.

    `act(p(X))` and `p(X)` are meant to describe identically, so that they
    intern to one species and become one database node: the projection has no
    edge between the two, and keeping them apart severs every upstream ->
    downstream path that runs through one and out of the other (6889 of them at
    two hops in the AD KG, which three hops cannot recover). Making them one
    node means the *subject* has to carry the active border too, so being active
    stops being a property of the term alone and becomes this set.

    Pass **every loaded term**, not just the projected roots: a protein is drawn
    active whenever BEL asserts any activity of it, including in a statement
    that mentions only the abundance, and including where the activity term
    itself is not projected.
    """
    subject_term_ids = set()
    for term in terms:
        if term.node_class != "Activity":
            continue
        for member in term.members():
            subject_term_ids.add(id(member))
    return frozenset(subject_term_ids)


def describe_term(term, cache=None, active_term_ids=frozenset()):
    """The :class:`TermDescription` of one BEL term. Pure, DB-free, memoised.

    `cache` is `{id(term): TermDescription}` and is shared across a whole build,
    so a term reached as a projected node and as a complex member is described
    once. `active_term_ids` is what
    :func:`collect_activity_subject_term_ids` returns: a term in it is described
    active, exactly as an `act(...)` of it is. Raises on a cycle, on a BEL class
    with no CellDesigner counterpart and on a modifier this module cannot route.
    """
    return _describe_term(
        term, {} if cache is None else cache, set(), active_term_ids
    )


def _describe_term(term, cache, in_progress, active_term_ids):
    if id(term) in cache:
        return cache[id(term)]
    if id(term) in in_progress:
        raise RuntimeError(
            f"BEL node {term.node_id} ({term.bel!r}) is a structural constituent of "
            "itself: the term graph has a cycle."
        )
    in_progress.add(id(term))
    try:
        description = _describe_term_uncached(
            term, cache, in_progress, active_term_ids
        )
    finally:
        in_progress.discard(id(term))
    cache[id(term)] = description
    return description


def _describe_term_uncached(term, cache, in_progress, active_term_ids):
    if term.node_class == "Activity":
        # `act(p(X))` *is* X in an active state: it takes its subject's class and
        # its subject's members, and the `ma()` code is dropped. Two activities
        # of one subject therefore describe identically and collapse into one
        # species -- accepted, and measured at 75 of the AD KG's 549. The subject
        # is described active too (`active_term_ids`), so `act(p(X))` and `p(X)`
        # collapse into *the same* species rather than two.
        members = term.members()
        if len(members) != 1:
            raise RuntimeError(
                f"BEL activity {term.node_id} ({term.bel!r}) has {len(members)} "
                "subject sub-terms; exactly one is expected."
            )
        if term.modifiers():
            raise RuntimeError(
                f"BEL activity {term.node_id} ({term.bel!r}) carries a modifier "
                "sub-term, which describe_term does not know how to place."
            )
        subject = _describe_term(members[0], cache, in_progress, active_term_ids)
        return dataclasses.replace(subject, active=True)

    species_class, template_class = _resolve_cd_classes(term)
    members = tuple(term.members())
    unknown = [
        edge_type
        for edge_type, _ in term.parts
        if edge_type not in MEMBER_EDGE_TYPES and edge_type not in MODIFIER_EDGE_TYPES
    ]
    if unknown:
        raise RuntimeError(
            f"BEL node {term.node_id} ({term.bel!r}) has {sorted(set(unknown))} "
            "sub-term edges, which are not entity members nor modifiers. Only BEL "
            "process terms carry those, and none is projected."
        )

    modifications = []
    structural_states = []
    name_suffix_parts = []
    compartment_name = None
    for edge_type, part in term.modifiers():
        if edge_type == "HAS__LOCATION":
            # `compartment` lives on the `Species` base, so a location needs no
            # per-class fallback -- it is the one modifier that always applies.
            location_name = part.name or part.bel
            if compartment_name is not None and location_name != compartment_name:
                raise RuntimeError(
                    f"BEL node {term.node_id} ({term.bel!r}) has more than one "
                    "location; a species has a single compartment."
                )
            compartment_name = location_name
            continue
        kind, payload = _describe_modifier(edge_type, part)
        if kind == "fragment" and issubclass(species_class, PROTEIN_SPECIES_CLASSES):
            # The class change happens even for a `frag("?")`: *that* is what says
            # the species is a fragment, and it says it whether or not the range
            # is known -- so the `None` payload below adds nothing on top.
            species_class, template_class = _TRUNCATED_CLASSES
            if payload is not None:
                structural_states.append(payload)
        elif payload is None:
            # The sub-term carries nothing to draw; never invent a placeholder.
            continue
        elif kind == "modification" and issubclass(
            species_class, PROTEIN_SPECIES_CLASSES
        ):
            modifications.append(payload)
        elif kind == "structural_state" and issubclass(
            species_class, PROTEIN_SPECIES_CLASSES + COMPLEX_SPECIES_CLASSES
        ):
            structural_states.append(payload)
        else:
            # A gene, an RNA, an abundance or a phenotype: no badge container, so
            # the modifier is folded into the species *name* -- never into the
            # template name, or one protein would declare N `<protein>` entries.
            name_suffix_parts.append(_name_suffix_part(kind, payload))

    if term.node_class == "Composite":
        # CellDesigner has no composite; the state says what the complex stands
        # for so a reader is not misled into thinking the members are bound.
        structural_states.append("composite")

    identifier = term.name
    if not identifier:
        if members:
            identifier = ":".join(
                sorted(
                    _describe_term(
                        member, cache, in_progress, active_term_ids
                    ).species_name
                    for member in members
                )
            )
        else:
            identifier = term.bel or term.node_id
    species_name = identifier
    if name_suffix_parts:
        species_name = f"{identifier} [{'|'.join(sorted(set(name_suffix_parts)))}]"
    return TermDescription(
        species_class=species_class,
        template_class=template_class,
        template_name=identifier if template_class is not None else None,
        species_name=species_name,
        active=id(term) in active_term_ids,
        compartment_name=compartment_name,
        # Sorted through `_residue_sort_key`: a residue name may be `None` (a
        # `pmod(Ph)` naming no site), and `None` does not compare with a string.
        modifications=tuple(
            sorted(
                set(modifications),
                key=lambda modification: (
                    _residue_sort_key(modification[0]),
                    modification[1].value,
                ),
            )
        ),
        structural_states=tuple(sorted(set(structural_states))),
        members=members,
    )


def _resolve_cd_classes(term):
    if term.node_class == "Abundance":
        if term.namespace in ABUNDANCE_MOLECULE_NAMESPACES:
            return momapy.celldesigner.SimpleMolecule, None
        return momapy.celldesigner.Unknown, None
    classes = BEL_CLASS_TO_CD_CLASSES.get(term.node_class)
    if classes is None:
        raise RuntimeError(
            f"BEL node {term.node_id} ({term.bel!r}) has class {term.node_class!r}, "
            "which has no CellDesigner counterpart. Add it to "
            "BEL_CLASS_TO_CD_CLASSES."
        )
    return classes


# ---------------------------------------------------------------------------
# templates, residues and compartments: the interning pass
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class BuildContext:
    """Everything :func:`make_species` needs, built once per run.

    `templates` and `residues` carry invariants (a) and (b), `compartments` and
    `compartment_roots` invariant (d); see the module docstring. `descriptions`
    is the shared `describe_term` memo, so it covers the projected roots *and*
    every term reached as a member.
    """

    descriptions: dict  # id(BelTerm) -> TermDescription
    templates: dict  # (template class, name) -> template object
    residues: dict  # (template class, name) -> {residue name: residue object}
    compartments: dict  # (collection name, location name) -> Compartment
    compartment_roots: dict  # collection name -> Compartment
    stats: dict


def _collect_terms(terms):
    """`terms` plus every term reachable through their members, depth first."""
    collected = {}
    stack = list(terms)
    while stack:
        term = stack.pop()
        if id(term) in collected:
            continue
        collected[id(term)] = term
        stack.extend(term.members())
    return list(collected.values())


def make_build_context(terms, active_term_ids=frozenset()):
    """The :class:`BuildContext` for the given root terms. Three passes, in order.

    A species is frozen, so its template cannot be patched once built: every
    residue a modification can name has to be known before the first template
    exists. Hence

    1. describe every term that will be materialised -- the roots **and** their
       recursive members;
    2. accumulate `{(template class, name): {residue name}}` over those
       descriptions, build one residue object per `(name, order)` and one
       template per `(class, name)`, and intern the template by value;
    3. **re-derive** the residue lookup from the *interned* template -- invariant
       (b).

    The interning by `(class, name)` is what keeps **one template object per
    protein within the export**, which invariant (b) depends on.

    Identity with the *stored* collections' templates is not this function's job:
    the export is standalone, and hash integration with a seeded
    `object_key_to_node` makes two value-equal elements one database node across
    save calls (see `2_10`).

    `active_term_ids` is threaded into every description; see
    :func:`collect_activity_subject_term_ids`.
    """
    all_terms = _collect_terms(terms)
    descriptions = {}
    for term in all_terms:
        _describe_term(term, descriptions, set(), active_term_ids)

    residue_names = {}
    for term in all_terms:
        description = descriptions[id(term)]
        if description.template_class is None:
            continue
        key = (description.template_class, description.template_name)
        names = residue_names.setdefault(key, set())
        if description.template_class in TEMPLATE_RESIDUE_FIELDS:
            names.update(residue_name for residue_name, _ in description.modifications)

    templates = {}
    residues = {}
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
        # One object per key, and the key is `(class, name)`: that is invariant
        # (a) within the export.
        template = template_class(
            name=template_name, **{field_name: template_residues}
        )
        templates[key] = template
        # Invariant (b): from the object the species will actually carry.
        residues[key] = {
            residue.name: residue for residue in getattr(template, field_name)
        }

    compartment_roots = {}
    compartments = {}
    for term in all_terms:
        description = descriptions[id(term)]
        if description.compartment_name is None:
            continue
        root = compartment_roots.get(term.collection_name)
        if root is None:
            # Undrawn, like the stored maps' `default` root: it exists to give
            # every BEL compartment an `outside`, which is what makes
            # `pd2af.utils._build_dot_graph` attach their clusters at all, and
            # what keeps a BEL `nucleus` value-distinct from a stored one.
            root = momapy.celldesigner.Compartment(name=term.collection_name)
            compartment_roots[term.collection_name] = root
        compartment_key = (term.collection_name, description.compartment_name)
        if compartment_key not in compartments:
            compartments[compartment_key] = momapy.celldesigner.Compartment(
                name=description.compartment_name, outside=root
            )

    stats = {
        "n_terms_described": len(all_terms),
        "n_templates": len(templates),
        "n_compartments": len(compartments),
        "n_compartment_roots": len(compartment_roots),
    }
    return BuildContext(
        descriptions=descriptions,
        templates=templates,
        residues=residues,
        compartments=compartments,
        compartment_roots=compartment_roots,
        stats=stats,
    )


# ---------------------------------------------------------------------------
# the species
# ---------------------------------------------------------------------------


def make_species(term, context, is_subunit=False, record=None):
    """One frozen momapy species for one BEL term. Always a **fresh object tree**.

    Invariant (c): the writer indexes complex subunits by `id()` and skips any
    top-level species whose `id()` is in that index, so memoising a species by
    node id -- the obvious way to get value-collapse -- would make the own glyph
    of every protein that is also a complex member disappear, and reroute every
    modulation targeting it to the enclosing complex. Interning happens only at
    the top level, on the finished object, in `commute_dm.bel_export`.

    A subunit keeps no compartment: CellDesigner puts an included species inside
    its complex, not inside a compartment box, and only top-level species
    contribute to `model.compartments`.

    When `record` is a list, `(term.node_id, species)` is appended for **every**
    species built, subunits included. That is what lets a subunit be annotated:
    a subunit is a per-occurrence object which `species_by_node_id` cannot
    recover, and 16 of the interface's UniProt identifiers are carried only by
    complex members.
    """
    description = context.descriptions[id(term)]
    fields = {
        "name": description.species_name,
        "active": description.active,
    }
    if not is_subunit and description.compartment_name is not None:
        fields["compartment"] = context.compartments[
            (term.collection_name, description.compartment_name)
        ]
    if description.template_class is not None:
        key = (description.template_class, description.template_name)
        fields["template"] = context.templates[key]
        residues = context.residues[key]
        if description.template_class in TEMPLATE_RESIDUE_FIELDS:
            fields["modifications"] = frozenset(
                momapy.celldesigner.Modification(
                    residue=residues[residue_name], state=state
                )
                for residue_name, state in description.modifications
            )
        else:
            fields["modifications"] = frozenset()
    if issubclass(
        description.species_class,
        PROTEIN_SPECIES_CLASSES + COMPLEX_SPECIES_CLASSES,
    ):
        fields["structural_states"] = frozenset(
            momapy.celldesigner.StructuralState(value=value)
            for value in description.structural_states
        )
    if issubclass(description.species_class, COMPLEX_SPECIES_CLASSES):
        fields["subunits"] = frozenset(
            make_species(member, context, is_subunit=True, record=record)
            for member in description.members
        )
    species = description.species_class(**fields)
    if record is not None:
        record.append((term.node_id, species))
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
