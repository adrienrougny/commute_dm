"""Comorbidity sub-maps out of the **stored** activity-flow CellDesigner maps.

The DB holds activity-flow maps (`COVID_DM_CD_AF`, `PD_DM_CD_AF`) produced by
`pd2af` in `keep-reactions` mode. An AF map has 0 reactions and its
`Modulation`s point straight at `Species`, so it already *contains*, as stored
objects, everything a sub-map needs: species, signed modulations, boolean logic
gates, their glyphs and arcs, and the `LayoutModelMapping` tying them together.

Four steps, of which only the third knows where a source came from:

1. the interface (`commute_dm.core.get_interface`) gives the seeds;
2. :func:`load_signed_influences` pulls the influence structure out of the
   database as an **edge list of node ids** and :class:`Influences` walks it in
   Python -- the walk knows nothing about CellDesigner or BEL;
3. :func:`load_collections_as_map` materialises the drawing material, once per
   run, as one merged `CellDesignerMap`;
4. :func:`make_submap_from_model_elements` assembles one sub-map out of selected
   model elements, reusing pd2af's layout helpers.

**Why this is so much smaller than the `ig.py` it replaces.** Collections are
imported with `integration_mode="hash"`, so two content-equal elements are *one*
database node; hydrating everything through one shared `node_id_to_object` cache
therefore makes that one node exactly one Python object. Object identity --
which `LayoutModelMapping` is keyed on, and which the CellDesigner writer
resolves aliases by -- then holds for free, across submap boundaries and across
collections. That single fact is what removes the old interning, canonical
element/glyph tables and gate-input re-pointing: there is nothing left to
canonicalise.

Geometry is *not* reused: stored arc positions mean nothing once maps are
merged, so callers must run `pd2af.utils.make_auto_layout` on the result. It
recomputes every position and segment, and fits the root layout.
"""

import collections
import dataclasses
import re

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.geometry
import pd2af.celldesigner.building_layout
import pd2af.celldesigner.building_model
import pd2af.utils

import commute_dm.queries


# ---------------------------------------------------------------------------
# the influence structure, as node ids
# ---------------------------------------------------------------------------

# The modulation classes whose sign is known. Matching `:Modulation` already
# leaves out every `Unknown*` class -- they are siblings of `Modulation` under
# `KnownOrUnknownModulation`, not subclasses -- and this leaves out the plain
# `Modulation`s, which carry no sign. `Inhibition` belongs here even though
# pd2af never builds one: the reader's `normalize_modulation_class` turns a
# negative modulation targeting a `Phenotype` back into an `Inhibition`, so
# reading pd2af's own output gives them. `Catalysis` and `PhysicalStimulation`
# cannot occur in an activity-flow map but do in the process-description ones,
# and the collections are a parameter. pylpg labels a node with the name of its
# class and of every ancestor, so the same tuple serves the Cypher label test
# and the `isinstance` test in `make_submap_from_model_elements`.
SIGNED_MODULATION_CLASSES = (
    momapy.celldesigner.PositiveInfluence,
    momapy.celldesigner.NegativeInfluence,
    momapy.celldesigner.Triggering,
    momapy.celldesigner.Inhibition,
    momapy.celldesigner.PhysicalStimulation,
    momapy.celldesigner.Catalysis,
)

_SIGNED_INFLUENCES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
    -[:HAS_MODEL]->(:CellDesignerModel)-[:HAS_MODEL_ELEMENT]->(modulation:Modulation)
WHERE collection.name IN $collection_names
  AND any(label IN labels(modulation) WHERE label IN $signed_classes)
MATCH (modulation)-[:HAS_SOURCE]->(source), (modulation)-[:HAS_TARGET]->(target)
RETURN DISTINCT elementId(source) AS source_node_id, elementId(target) AS target_node_id
"""

# The modulations the query above leaves out, by concrete class. Not used by the
# walk: it is the guard that replaces `ig._classify_modulation`, which *raised*
# on an unrecognised class rather than skipping it -- silently skipping is how
# `make_ig_in_db` rotted into dead code after momapy renamed things under it.
_EXCLUDED_INFLUENCES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
    -[:HAS_MODEL]->(:CellDesignerModel)-[:HAS_MODEL_ELEMENT]->(modulation:KnownOrUnknownModulation)
WHERE collection.name IN $collection_names
  AND NOT any(label IN labels(modulation) WHERE label IN $signed_classes)
RETURN DISTINCT elementId(modulation) AS modulation_node_id,
       labels(modulation) AS modulation_labels
"""

# The species of the collections, so that a selection can be narrowed to the
# nodes that can carry an annotation. `HAS_MODEL_ELEMENT` is deliberate: it
# reaches the 1344 complex-subunit species as well as the 3540 top-level ones,
# and a subunit does carry annotations.
_SPECIES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
    -[:HAS_MODEL]->(:CellDesignerModel)-[:HAS_MODEL_ELEMENT]->(species:Species)
WHERE collection.name IN $collection_names
RETURN collect(DISTINCT elementId(species)) AS species_node_ids
"""

# A gate's inputs go through a `SimpleSpeciesReference`, not straight to the
# species -- getting this path wrong is what made the old gate queries match 0
# rows.
_GATE_INPUTS_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
    -[:HAS_MODEL]->(:CellDesignerModel)-[:HAS_MODEL_ELEMENT]->(gate:BooleanLogicGate)
WHERE collection.name IN $collection_names
MATCH (gate)-[:HAS_INPUT]->(:SimpleSpeciesReference)-[:HAS_REFERRED_ELEMENT]->(species:Species)
RETURN elementId(gate) AS gate_node_id,
       collect(DISTINCT elementId(species)) AS input_node_ids
"""

# Concrete classes the sign filter is *expected* to leave out. Anything else
# raises: a new modulation class must be classified deliberately, here or in
# `SIGNED_MODULATION_CLASSES`.
KNOWN_UNSIGNED_MODULATION_CLASSES = frozenset(
    {
        "Modulation",
        "UnknownModulation",
        "UnknownPositiveInfluence",
        "UnknownNegativeInfluence",
        "UnknownInhibition",
        "UnknownTriggering",
        "UnknownPhysicalStimulation",
        "UnknownCatalysis",
    }
)

# Labels every modulation carries whatever its class, so they say nothing about
# it and are ignored when checking that an excluded class was excluded on
# purpose.
_STRUCTURAL_LABELS = frozenset(
    {
        "BaseNode",
        "MapElement",
        "ModelElement",
        "CellDesignerModelElement",
        "KnownOrUnknownModulation",
        "Modulation",
    }
)


@dataclasses.dataclass(frozen=True)
class Influences:
    """The signed influence structure of one or more collections, as node ids.

    Nodes are species **and** boolean logic gates, so an influence through a
    gate is two steps: the gate's inputs influence the gate, and the gate
    influences whatever its modulations target.

    Holds no momapy objects. `node_id` throughout is the Neo4j `elementId`,
    which is also the key of momapy_kb's `node_id_to_object` hydration cache --
    that is what turns a selection back into model elements.

    One `Influences` covers both collections, as one `AfIndex` used to. Loading
    one per collection instead would be equivalent, measurably so: the two AF
    influence graphs share **no** node, so an upstream walk from a COVID seed
    cannot wander into PD anyway. (201 species nodes *are* in both collections,
    but all 201 are complex subunits and none is a modulation endpoint or a gate
    input.) Keeping one object is simply less to thread through.
    """

    influencing_node_ids: dict
    influenced_node_ids: dict
    gate_input_node_ids: dict
    species_node_ids: frozenset

    def upstream(self, seed_node_ids, max_level):
        """The elements that influence the seeds, at most `max_level` steps away."""
        return _select_within_levels(
            self.influencing_node_ids, seed_node_ids, max_level
        )

    def downstream(self, seed_node_ids, max_level):
        """The elements the seeds influence, at most `max_level` steps away."""
        return _select_within_levels(self.influenced_node_ids, seed_node_ids, max_level)

    def close_over_gates(self, node_ids):
        """Add the inputs of every selected boolean logic gate.

        A gate whose inputs are not in the map refers to species the map does
        not contain, which the writer emits as a dangling input.
        """
        node_ids = set(node_ids)
        return node_ids.union(
            *(
                self.gate_input_node_ids[node_id]
                for node_id in node_ids & self.gate_input_node_ids.keys()
            )
        )

    def species_only(self, node_ids):
        """The species among `node_ids`.

        Drops the boolean logic gates: a gate carries no annotation, so it
        cannot contribute to a gene set.
        """
        return set(node_ids) & self.species_node_ids


def load_signed_influences(session, collection_names):
    """The signed influence graph of the given collections. See :class:`Influences`.

    Raises when a modulation is neither signed nor a known unsigned class, so
    that a momapy addition is a loud failure rather than a silently smaller
    graph.
    """
    parameters = {
        "collection_names": list(collection_names),
        "signed_classes": [class_.__name__ for class_ in SIGNED_MODULATION_CLASSES],
    }
    for row in session.execute_query(_EXCLUDED_INFLUENCES_QUERY, parameters):
        concrete = [
            label
            for label in row["modulation_labels"]
            if label not in _STRUCTURAL_LABELS
        ]
        unexpected = [
            label for label in concrete if label not in KNOWN_UNSIGNED_MODULATION_CLASSES
        ]
        if unexpected:
            raise RuntimeError(
                f"modulation {row['modulation_node_id']} has unrecognised class "
                f"{sorted(unexpected)}. Add it to SIGNED_MODULATION_CLASSES or to "
                "KNOWN_UNSIGNED_MODULATION_CLASSES."
            )
    species_node_ids = frozenset(
        session.execute_query(_SPECIES_QUERY, parameters)[0]["species_node_ids"]
    )
    gate_input_node_ids = {
        row["gate_node_id"]: set(row["input_node_ids"])
        for row in session.execute_query(_GATE_INPUTS_QUERY, parameters)
    }
    influencing_node_ids = collections.defaultdict(set)
    influenced_node_ids = collections.defaultdict(set)
    for row in session.execute_query(_SIGNED_INFLUENCES_QUERY, parameters):
        influencing_node_ids[row["target_node_id"]].add(row["source_node_id"])
        influenced_node_ids[row["source_node_id"]].add(row["target_node_id"])
    for gate_node_id, input_node_ids in gate_input_node_ids.items():
        influencing_node_ids[gate_node_id] |= input_node_ids
        for input_node_id in input_node_ids:
            influenced_node_ids[input_node_id].add(gate_node_id)
    return Influences(
        influencing_node_ids=dict(influencing_node_ids),
        influenced_node_ids=dict(influenced_node_ids),
        gate_input_node_ids=gate_input_node_ids,
        species_node_ids=species_node_ids,
    )


def _select_within_levels(neighbour_node_ids, seed_node_ids, max_level):
    """Breadth-first search bounded by `max_level`; a negative bound means none.

    The seeds are always selected. `max_level` means hops, plainly.
    """
    selected = set(seed_node_ids)
    frontier = set(seed_node_ids)
    level = 0
    while frontier and (max_level < 0 or level < max_level):
        frontier = {
            neighbour_node_id
            for node_id in frontier
            for neighbour_node_id in neighbour_node_ids.get(node_id, ())
        } - selected
        selected |= frontier
        level += 1
    return selected


# ---------------------------------------------------------------------------
# the drawing material
# ---------------------------------------------------------------------------

_COLLECTION_MAPS_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(map:CellDesignerMap)
WHERE collection.name IN $collection_names
RETURN DISTINCT map
"""


def load_collections_as_map(session, collection_names, node_id_to_object):
    """Every map of the given collections, merged into one (a few minutes).

    `node_id_to_object` is momapy_kb's hydration cache. It is what turns the
    node ids the walk returns into model elements, and sharing it across the
    whole run is what makes the merge sound -- see the module docstring.

    `prewarm_session` is a momapy_kb requirement, not something this module
    introduces: `execute_query_as_objects` can only rebuild an object whose
    momapy class has a registered pylpg node class, and the session registers
    only what it saved or what a static type hint mentions -- concrete
    subclasses like `GenericProteinLayout` are missed.
    """
    commute_dm.queries.prewarm_session(session)
    return merge_maps(
        [
            row[0]
            for row in session.execute_query_as_objects(
                _COLLECTION_MAPS_QUERY,
                {"collection_names": list(collection_names)},
                node_id_to_object=node_id_to_object,
            )
        ]
    )


def merge_maps(cd_maps):
    """One CellDesignerMap holding the model elements, the layout elements and
    the mapping entries of all the given maps.

    A species drawn in two submaps is one node in the database and therefore one
    `Species` here, which is what lets a walk cross a submap boundary and what
    makes merging the mappings lossless.
    """

    def union(attribute):
        return frozenset().union(
            *(getattr(cd_map.model, attribute) for cd_map in cd_maps)
        )

    mapping_items = {}
    for cd_map in cd_maps:
        if cd_map.layout_model_mapping is None:
            continue
        mapping_items.update(cd_map.layout_model_mapping)
    return momapy.celldesigner.CellDesignerMap(
        model=momapy.celldesigner.CellDesignerModel(
            compartments=union("compartments"),
            species=union("species"),
            species_templates=union("species_templates"),
            boolean_logic_gates=union("boolean_logic_gates"),
            modulations=union("modulations"),
        ),
        # This map is a lookup structure and is never rendered, but
        # `Node.position`, `Node.width` and `Node.height` have no defaults, so
        # they have to be given.
        layout=momapy.celldesigner.CellDesignerLayout(
            position=momapy.geometry.Point(0.0, 0.0),
            width=0.0,
            height=0.0,
            layout_elements=tuple(
                layout_element
                for cd_map in cd_maps
                for layout_element in cd_map.layout.layout_elements
            ),
        ),
        layout_model_mapping=momapy.core.mapping.LayoutModelMapping(mapping_items),
    )


def count_mapping_entries(cd_maps):
    """`(merged, sum of the parts)` mapping entry counts, for verification.

    A `LayoutModelMapping`'s forward dict is equality-keyed, so two
    content-equal layout keys coming from two maps would collapse and one
    glyph -> model element entry would be lost. Hash integration makes them the
    same object and the counts equal; this is the check that says so.
    """
    merged = merge_maps(cd_maps)
    return len(merged.layout_model_mapping), sum(
        len(cd_map.layout_model_mapping or {}) for cd_map in cd_maps
    )


# ---------------------------------------------------------------------------
# identity invariants
# ---------------------------------------------------------------------------
#
# momapy elements are equal by value with `id_` excluded, and both this module
# and the CellDesigner writer rely on *object identity* in places where a
# value-equal duplicate would go unnoticed until the written file is read back.
# The four checks below are run on every sub-map just before `renumber_ids`, so
# a violation is a loud failure here rather than a `KeyError` in the reader.
#
# They are cheap (a walk of the sub-map, which holds tens of elements) and they
# hold for stored-only data too -- hash integration makes one database node
# exactly one Python object -- so they are not BEL-specific; the BEL side is
# simply the one that has to *establish* them (see `commute_dm.bel2cd`).


def iter_subunits(species):
    """Every subunit of `species`, at any depth."""
    for subunit in getattr(species, "subunits", ()) or ():
        yield subunit
        yield from iter_subunits(subunit)


def iter_species_and_subunits(species_list):
    """Every species of `species_list` and, recursively, its subunits.

    A complex's members are not top-level species, so `model.species` alone
    misses them -- and they are where half the modifications and templates of a
    BEL sub-map live.
    """
    for species in species_list:
        yield species
        yield from iter_subunits(species)


def _get_template_residues(template):
    """A protein template's `modification_residues`, a gene/RNA one's `regions`."""
    for field_name in ("modification_residues", "regions"):
        container = getattr(template, field_name, None)
        if container is not None:
            return container
    return ()


def check_species_template_identity(cd_map):
    """One template *object* per template value. Invariant (a).

    `make_submap_from_model_elements` derives `model.species_templates` by value
    via `collect_templates_from_species`, so two value-equal but distinct
    template objects collapse into one member of that frozenset while each
    species keeps pointing at its own object. The writer then emits a single
    `<protein>` and the other species emit a `<proteinReference>` to an id that
    was never written -- the file is written without error and fails to read back
    with a `KeyError`.
    """
    identities = collections.defaultdict(set)
    for species in iter_species_and_subunits(cd_map.model.species):
        template = getattr(species, "template", None)
        if template is not None:
            identities[template].add(id(template))
    duplicated = [
        template for template, ids in identities.items() if len(ids) > 1
    ]
    if duplicated:
        raise RuntimeError(
            f"{len(duplicated)} species template(s) exist as several distinct "
            f"objects with the same value, e.g. {duplicated[0].name!r}. They would "
            "collapse into one `<protein>` declaration while the other species "
            "reference an id that is never written. Intern templates by value."
        )


def check_modification_residue_identity(cd_map):
    """A modification's residue *is* an object of its template's container. Invariant (b).

    The writer emits `<modification residue=...>` from the modification and
    `<modificationResidue id=...>` from the template independently; after
    `renumber_ids` two value-equal residue objects carry two different ids, so
    the modification names a residue the list does not declare and the reader's
    unguarded lookup raises.
    """
    for species in iter_species_and_subunits(cd_map.model.species):
        template = getattr(species, "template", None)
        modifications = getattr(species, "modifications", ()) or ()
        if not modifications:
            continue
        residue_ids = {id(residue) for residue in _get_template_residues(template)}
        for modification in modifications:
            residue = modification.residue
            if residue is None:
                raise RuntimeError(
                    f"a modification of species {species.name!r} has no residue; the "
                    "reader looks a modification's residue up unguarded."
                )
            if id(residue) not in residue_ids:
                raise RuntimeError(
                    f"species {species.name!r} has a modification on residue "
                    f"{residue.name!r} that is not an object of its template's "
                    "residue container. The written modification would name a "
                    "residue that is never declared."
                )


def find_subunit_top_level_aliases(cd_map):
    """`[(complex, subunit)]` where a subunit object *is* a top-level species.

    The writer keys `build_subunit_to_complex` by `id(subunit)` and skips any
    species whose `id()` is in that index, so such a species is written as
    neither `<species>` nor `<speciesAlias>` while its modulations still point at
    it -- a file that writes cleanly and fails to read back with `KeyError`.

    This is **not** a BEL problem. Hash integration makes one database node
    exactly one Python object, so a stored species that the walk selects *and*
    that a selected complex holds as a subunit is one object playing both roles;
    it reaches a selection as a boolean-gate input or as an interface seed.
    Measured on COVID x AD before the fix: 7 such pairs over 5 of 101 maps
    (`CASP1` in `NLRP3 oligomer:PYCARD:CASP1`, `TRAF2` in `TRAF2:ERN1:unfolded
    protein`), and those 5 were exactly the 5 that failed to read back. It is
    latent in COVID -> PD too; that pairing's maps simply never select such a
    pair.

    :func:`make_submap_from_model_elements` resolves it by giving the standalone
    occurrence its own object, so this must find nothing -- which is what
    :func:`check_identity_invariants` asserts.
    """
    top_level_ids = {id(species) for species in cd_map.model.species}
    aliases = []
    for complex_ in cd_map.model.species:
        for subunit in iter_subunits(complex_):
            if id(subunit) in top_level_ids:
                aliases.append((complex_, subunit))
    return aliases


def check_subunit_layouts(cd_map):
    """Every subunit is drawn inside its complex, and is its own object. Invariant (c).

    The writer emits an included species' alias only when
    `get_child_layout_elements(subunit, complex)` returns a layout that is also
    in the complex layout's `layout_elements`, so both the nesting and the
    mapping entry are needed. And no subunit object may *be* a top-level species
    -- see :func:`find_subunit_top_level_aliases`.
    """
    aliases = find_subunit_top_level_aliases(cd_map)
    if aliases:
        complex_, subunit = aliases[0]
        raise RuntimeError(
            f"{len(aliases)} subunit(s) are the same object as a top-level species, "
            f"e.g. {subunit.name!r} of complex {complex_.name!r}. The writer would "
            "emit neither a species nor an alias for the standalone occurrence, and "
            "the file would fail to read back."
        )
    for complex_ in cd_map.model.species:
        subunits = getattr(complex_, "subunits", ()) or ()
        if not subunits:
            continue
        complex_layout_elements = [
            key
            for key in cd_map.get_mapping(complex_) or []
            if not isinstance(key, frozenset)
        ]
        for subunit in subunits:
            drawn = [
                layout_element
                for complex_layout in complex_layout_elements
                for layout_element in cd_map.layout_model_mapping
                .get_child_layout_elements(subunit, complex_)
                if layout_element in complex_layout.layout_elements
            ]
            if not drawn:
                raise RuntimeError(
                    f"subunit {subunit.name!r} of complex {complex_.name!r} has no "
                    "mapped layout element nested in the complex's glyph, so the "
                    "writer emits no alias for it."
                )


def check_compartment_identity(cd_map):
    """One compartment *object* per compartment value. Invariant (d).

    The same collapse as invariant (a), one level down: the writer emits one
    `<compartment>` per member of the model's compartment frozenset, so two
    value-equal compartment objects leave a species referencing an id that was
    never written.
    """
    identities = collections.defaultdict(set)
    for species in cd_map.model.species:
        compartment = species.compartment
        if compartment is not None:
            identities[compartment].add(id(compartment))
    for compartment in cd_map.model.compartments:
        identities[compartment].add(id(compartment))
        if compartment.outside is not None:
            identities[compartment.outside].add(id(compartment.outside))
    duplicated = [
        compartment for compartment, ids in identities.items() if len(ids) > 1
    ]
    if duplicated:
        raise RuntimeError(
            f"{len(duplicated)} compartment(s) exist as several distinct objects "
            f"with the same value, e.g. {duplicated[0].name!r}. One of them would "
            "not be written, and the species inside it would reference an unknown "
            "compartment id."
        )


def check_identity_invariants(cd_map):
    """The four identity invariants, in one call."""
    check_species_template_identity(cd_map)
    check_modification_residue_identity(cd_map)
    check_subunit_layouts(cd_map)
    check_compartment_identity(cd_map)


# ---------------------------------------------------------------------------
# label measuring and fitted synthetic glyphs
# ---------------------------------------------------------------------------

# These live here rather than in `commute_dm.bel2cd` because `commute_dm.core`
# needs `make_fitted_synthetic_layout` for every sub-map's synthetic central
# node, in both pairings, while `bel2cd` is export-only. `bel2cd` imports
# them back; this module imports nothing of it, so the dependency is acyclic.

# Wrapping keeps a glyph from becoming absurdly wide. Bare HGNC symbols are
# short, but a CHEBI abundance name is not -- the longest in the AD KG runs to
# 130 characters -- and a complex is named after its members. 260pt holds a
# whole `beta-D-GalNAc-(1->4)-` chunk on one line.
_LABEL_MAX_WIDTH = 260.0
# Breathing room between the text and the glyph outline, on both axes.
_LABEL_PADDING = 12.0
# Long labels break after a `,`, which separates a complex's members, so a break
# there leaves whole names on a line. The trailing alternative keeps the final
# chunk.
_LABEL_CHUNK_PATTERN = re.compile(r"[^,]*,|[^,]+")
_MEASURING_POSITION = momapy.geometry.Point(0.0, 0.0)


def _measure(text):
    """`(width, height)` of `text` as a default `TextLayout` would draw it."""
    bounding_box = momapy.core.layout.TextLayout(
        text=text, position=_MEASURING_POSITION
    ).bbox()
    return bounding_box.width, bounding_box.height


def _wrap_label_text(text, max_width=_LABEL_MAX_WIDTH):
    """`(wrapped_text, width, height)` for a label wrapped to `max_width`.

    `TextLayout` does **not** wrap on its own -- setting its `width` leaves the
    bounding box unchanged -- but it does honour newlines, so the wrapping is
    done here. The width budget is derived from one real measurement of the
    whole string rather than from a guessed per-character width, and the
    assembled result is measured again, so the returned size is the true one
    even where a break landed badly.
    """
    width, height = _measure(text)
    if width <= max_width:
        return text, width, height
    # `text` is never empty here: it measured wider than `max_width`.
    character_budget = max(int(len(text) * max_width / width), 8)
    chunks = []
    for chunk in _LABEL_CHUNK_PATTERN.findall(text):
        # A single chunk with no break opportunity (a long chemical name) still
        # has to be split, or the glyph would be as wide as the unwrapped label.
        while len(chunk) > character_budget:
            chunks.append(chunk[:character_budget])
            chunk = chunk[character_budget:]
        if chunk:
            chunks.append(chunk)
    lines = []
    current_line = ""
    for chunk in chunks:
        if current_line and len(current_line) + len(chunk) > character_budget:
            lines.append(current_line)
            current_line = chunk
        else:
            current_line += chunk
    if current_line:
        lines.append(current_line)
    wrapped_text = "\n".join(lines)
    width, height = _measure(wrapped_text)
    return wrapped_text, width, height


def make_fitted_synthetic_layout(species, index):
    """A synthetic glyph for `species`, sized to hold its (wrapped) label.

    `pd2af`'s `make_synthetic_layout` gives the right layout class, a throwaway
    position and a label, but leaves the class's default width and height.
    `make_auto_layout` *preserves* the size it is given -- its `_build_dot_graph`
    sets each dot node's size from the layout element, and only compartments are
    resized afterwards -- so sizing here is what reaches the output.

    Used for `commute_dm.core`'s synthetic central node, whose `Unknown` default
    of 60x30 is too small for some display names. BEL species go through
    `commute_dm.bel2cd.make_species_layouts`, which does the same fitting and
    then makes room for subunits and badges.
    """
    layout_element = pd2af.celldesigner.building_layout.make_synthetic_layout(
        species, index
    )
    wrapped_text, width, height = _wrap_label_text(layout_element.label.text)
    return dataclasses.replace(
        layout_element,
        label=dataclasses.replace(layout_element.label, text=wrapped_text),
        width=max(layout_element.width, width + _LABEL_PADDING),
        height=max(layout_element.height, height + _LABEL_PADDING),
    )


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------


def get_layout_element_for_model_element(cd_map, model_element, preferred_ids=None):
    """The layout element the map draws `model_element` with.

    A species and a boolean logic gate have exactly one. A compartment drawn in
    several of the merged maps has one per map, and the root `default`
    compartment has none; the lowest id keeps the choice deterministic.

    `preferred_ids` is a set of `id()`s to choose from first. It exists for the
    species that are drawn **both** standalone and as a complex subunit: those
    have two glyphs, and passing the ids of `cd_map.layout.layout_elements`
    picks the standalone one rather than the glyph nested inside the complex.
    """
    layout_elements = [
        key
        for key in cd_map.get_mapping(model_element) or []
        if not isinstance(key, frozenset)
    ]
    if preferred_ids:
        preferred = [
            layout_element
            for layout_element in layout_elements
            if id(layout_element) in preferred_ids
        ]
        layout_elements = preferred or layout_elements
    return min(
        layout_elements,
        key=lambda layout_element: layout_element.id_ or "",
        default=None,
    )


def get_arc_for_modulation(cd_map, modulation):
    """The arc the map draws `modulation` with.

    A modulation is mapped to a frozenset key holding the arc together with the
    layout elements of its source and target; the arc is the only one of the
    three that is not a node. `None` for a modulation the map does not draw.
    """
    for key in cd_map.get_mapping(modulation) or []:
        if isinstance(key, frozenset):
            for layout_element in key:
                if isinstance(layout_element, momapy.core.layout.Arc):
                    return layout_element
    return None


# The signed classes missing from
# `pd2af.celldesigner.building_layout._MODULATION_CLASS_TO_LAYOUT_CLASS`, which
# only covers what pd2af builds. These three cannot be built by pd2af but can be
# read back out of a stored map, so without them a modulation with no stored arc
# would fail on a `KeyError` rather than on anything informative.
_EXTRA_MODULATION_ARC_CLASSES = {
    momapy.celldesigner.Inhibition: momapy.celldesigner.InhibitionLayout,
    momapy.celldesigner.Catalysis: momapy.celldesigner.CatalysisLayout,
    momapy.celldesigner.PhysicalStimulation: (
        momapy.celldesigner.PhysicalStimulationLayout
    ),
}


def _make_modulation_arc(modulation, source_layout_element, target_layout_element):
    """An arc for a modulation the source map does not draw.

    In practice only the synthesized influences take this path.
    """
    arc_class = _EXTRA_MODULATION_ARC_CLASSES.get(type(modulation))
    if arc_class is None:
        return pd2af.celldesigner.building_layout.make_modulation_arc(
            modulation, source_layout_element, target_layout_element
        )
    return arc_class(
        source=source_layout_element,
        target=target_layout_element,
        segments=tuple(
            pd2af.utils.make_arc_segments_from_source_and_target(
                source_layout_element, target_layout_element
            )
        ),
    )


def make_submap_from_model_elements(
    source_map,
    model_elements,
    fills=None,
    extra_species=(),
    extra_influences=(),
):
    """A CellDesignerMap made of model elements selected from `source_map`.

    Each keeps the layout element `source_map` draws it with, recoloured per
    `fills` (a `{model element: momapy.coloring colour}` dict). `extra_species`
    is a list of `(species, layout_element)` with no counterpart in `source_map`
    -- the species standing for the interface protein today, materialised BEL
    nodes later -- and `extra_influences` a list of `(source, target,
    modulation_class)` added to the induced ones.

    Everything here is a frozen momapy object: the stored elements are reused as
    they are, and `dataclasses.replace` copies the few that need a different
    fill or different endpoints, sharing their children. The one builder round
    trip is :func:`renumber_ids`.

    The layout still carries the source maps' geometry, which means nothing once
    maps are merged: the caller runs `pd2af.utils.make_auto_layout`, which also
    fits the root layout.
    """
    fills = fills or {}
    layout_elements = []
    mapping = momapy.core.mapping.LayoutModelMappingBuilder()
    layout_element_of = {}

    model_elements = list(model_elements)
    # A selected species that a selected complex also holds as a subunit is
    # **one object**: hash integration makes one database node exactly one
    # Python object, and such a species reaches the selection as a boolean-gate
    # input or as an interface seed. The writer indexes subunits by `id()` and
    # skips any species whose `id()` is in that index, so it would emit neither
    # `<species>` nor `<speciesAlias>` for the standalone occurrence while its
    # modulations still point at it -- a file that writes cleanly and fails to
    # read back with `KeyError` (measured: exactly the 5 of 101 COVID x AD maps
    # that hit this). The standalone occurrence therefore gets its **own
    # object**, and its own glyph.
    #
    # The copy is value-*equal* to the original, which is what keeps the rest of
    # this function unchanged: `layout_element_of` and `fills` are equality-keyed,
    # so either object looks the two up. Only `model.species` membership and the
    # modulations' endpoints care about identity, and both are handled below.
    subunit_ids = {
        id(subunit)
        for model_element in model_elements
        for subunit in iter_subunits(model_element)
    }
    promoted = {
        id(model_element): dataclasses.replace(model_element)
        for model_element in model_elements
        if id(model_element) in subunit_ids
    }
    # The stored glyphs that are drawn at the top level of a map -- as opposed to
    # nested inside a complex's glyph. A species drawn both ways has one of each.
    top_level_layout_element_ids = {
        id(layout_element) for layout_element in source_map.layout.layout_elements
    }

    def resolve(model_element):
        return promoted.get(id(model_element), model_element)

    for model_element in model_elements:
        layout_element = get_layout_element_for_model_element(
            source_map, model_element, preferred_ids=top_level_layout_element_ids
        )
        if layout_element is None:
            raise RuntimeError(
                f"no layout element for {type(model_element).__name__} "
                f"{model_element.id_!r}: every species and boolean logic gate of a "
                "stored map is expected to be drawn exactly once"
            )
        resolved_model_element = resolve(model_element)
        if (
            resolved_model_element is not model_element
            and id(layout_element) not in top_level_layout_element_ids
        ):
            # The species has no standalone glyph, only the one nested in the
            # complex. Reusing that object would put the same layout element both
            # at the top level and inside the complex, and the equality-keyed
            # mapping would then hold one of the two model elements only. A
            # nudged position is enough to make the copy value-distinct; every
            # position here is thrown away by `make_auto_layout` anyway.
            layout_element = dataclasses.replace(
                layout_element,
                position=momapy.geometry.Point(
                    layout_element.position.x + 0.5, layout_element.position.y
                ),
            )
        fill = fills.get(model_element)
        if fill is not None:
            layout_element = dataclasses.replace(layout_element, fill=fill)
        layout_element_of[resolved_model_element] = layout_element
        mapping.add_mapping(layout_element, resolved_model_element)
        # `dataclasses.replace` is shallow, so the descendants below are the
        # stored objects either way and resolve in the stored mapping. The
        # source map already knows what the sub-elements of a layout element
        # stand for -- a complex's subunits, a species' modifications; without
        # those entries the writer refers to aliases it never writes.
        for descendant in layout_element.descendants():
            sub_model_element = source_map.get_mapping(descendant)
            if sub_model_element is not None and not isinstance(
                sub_model_element, list
            ):
                mapping.add_mapping(descendant, sub_model_element)
    for species, layout_element in extra_species:
        layout_element_of[species] = layout_element
        mapping.add_mapping(layout_element, species)

    species = {
        model_element
        for model_element in layout_element_of
        if isinstance(model_element, momapy.celldesigner.Species)
    }
    gates = {
        model_element
        for model_element in layout_element_of
        if isinstance(model_element, momapy.celldesigner.BooleanLogicGate)
    }
    # `(modulation as the sub-map holds it, the stored modulation it came from)`.
    # The pair is needed because a modulation whose endpoint was promoted above
    # has to be rebuilt on the promoted object -- the writer resolves a
    # participant by identity -- while the *stored* one is what
    # `get_arc_for_modulation` looks the reusable arc up by (identity again, so a
    # rebuilt copy would silently fall back to a synthesized arc).
    #
    # Deduplicated, and induced first so a modulation the source map draws keeps
    # its stored arc. `model.modulations` is a frozenset while the loop below
    # appends one arc per list element, so a duplicate would give two arcs for
    # one modulation -- and `extra_influences` can hold duplicates whenever
    # `node_id_to_object` is many-to-one, which it is on the BEL side (several
    # BEL node ids can stand for one species).
    modulation_pairs = {}
    for stored_modulation in source_map.model.modulations:
        if not isinstance(stored_modulation, SIGNED_MODULATION_CLASSES):
            continue
        if (
            stored_modulation.source not in layout_element_of
            or stored_modulation.target not in layout_element_of
        ):
            continue
        modulation = stored_modulation
        if (
            id(stored_modulation.source) in promoted
            or id(stored_modulation.target) in promoted
        ):
            modulation = dataclasses.replace(
                stored_modulation,
                source=resolve(stored_modulation.source),
                target=resolve(stored_modulation.target),
            )
        modulation_pairs.setdefault(modulation, (modulation, stored_modulation))
    for source, target, modulation_class in extra_influences:
        modulation = modulation_class(source=resolve(source), target=resolve(target))
        modulation_pairs.setdefault(modulation, (modulation, None))
    modulation_pairs = list(modulation_pairs.values())
    modulations = [modulation for modulation, _ in modulation_pairs]

    compartments = pd2af.celldesigner.building_model.collect_ancestor_compartments(
        {
            model_element.compartment
            for model_element in species
            if model_element.compartment is not None
        }
    )
    # `default` is the one id `renumber_ids` leaves alone (the writer and the
    # reader special-case it), so two root compartments carrying it would be
    # written as two `<compartment id="default">` and the reader would silently
    # keep one -- losing a compartment, and resolving one root's species to the
    # other.
    if len([c for c in compartments if c.id_ == "default"]) > 1:
        raise RuntimeError(
            "several root compartments with id 'default': the source maps disagree "
            "on the root compartment's attributes, so they were not merged on "
            "import. Make them identical and re-run 2_00_save_collections."
        )
    # Compartments first: they render an opaque interior, so a later one would
    # hide the species drawn inside it (pd2af orders them the same way).
    for compartment in pd2af.celldesigner.building_model.compartments_outermost_first(
        compartments
    ):
        layout_element = get_layout_element_for_model_element(source_map, compartment)
        if layout_element is None:  # the root `default` compartment is not drawn
            continue
        layout_elements.append(layout_element)
        # `make_auto_layout` looks compartments up through the mapping to build
        # its dot clusters, so an unmapped compartment glyph would be a box that
        # nothing is placed inside.
        mapping.add_mapping(layout_element, compartment)
    layout_elements.extend(layout_element_of.values())

    for modulation, stored_modulation in modulation_pairs:
        source_layout_element = layout_element_of[modulation.source]
        target_layout_element = layout_element_of[modulation.target]
        arc = (
            None
            if stored_modulation is None
            else get_arc_for_modulation(source_map, stored_modulation)
        )
        if arc is not None:
            # reused like a node; only its geometry is worthless once maps are
            # merged, and `make_auto_layout` rebuilds every segment anyway
            arc = dataclasses.replace(
                arc,
                source=source_layout_element,
                target=target_layout_element,
                segments=tuple(
                    pd2af.utils.make_arc_segments_from_source_and_target(
                        source_layout_element, target_layout_element
                    )
                ),
            )
        else:
            arc = _make_modulation_arc(
                modulation, source_layout_element, target_layout_element
            )
        layout_elements.append(arc)
        pd2af.celldesigner.building_layout.add_modulation_mapping(
            mapping, arc, source_layout_element, target_layout_element, modulation
        )
    # The CellDesigner writer locates a gate's inputs by scanning for logic arcs
    # sourced at the gate glyph, so without these every gate is written
    # input-less. pd2af adds no mapping entry for a logic arc either.
    for gate in gates:
        for gate_input in gate.inputs:
            layout_elements.append(
                pd2af.celldesigner.building_layout.make_logic_arc(
                    layout_element_of[gate],
                    layout_element_of[gate_input.referred_element],
                )
            )

    cd_map = momapy.celldesigner.CellDesignerMap(
        model=momapy.celldesigner.CellDesignerModel(
            compartments=frozenset(compartments),
            species=frozenset(species),
            species_templates=frozenset(
                pd2af.celldesigner.building_model.collect_templates_from_species(species)
            ),
            boolean_logic_gates=frozenset(gates),
            modulations=frozenset(modulations),
        ),
        # Throwaway geometry: `make_auto_layout` calls `harmonize_root_layout`,
        # which fits the root around its children.
        layout=momapy.celldesigner.CellDesignerLayout(
            position=momapy.geometry.Point(0.0, 0.0),
            width=0.0,
            height=0.0,
            layout_elements=tuple(layout_elements),
        ),
        layout_model_mapping=mapping.build(),
    )
    check_identity_invariants(cd_map)
    return renumber_ids(cd_map)


def renumber_ids(cd_map):
    """Give every element of the map a fresh `id_`.

    An id is unique within a source file, and a submap merges dozens of them, so
    ids collide -- species ids, compartment ids and layout element ids alike.
    The writer gives one xml id to every element sharing an id and the reader
    indexes by it, so two elements sharing an id are written and read back as
    one.

    The literal `default` compartment id is excluded: the writer and the reader
    special-case it, and renumbering it produces `KeyError: ''`.
    """
    object_to_builder = {}
    map_builder = momapy.builder.builder_from_object(
        cd_map, object_to_builder=object_to_builder
    )
    for n, element in enumerate(
        cd_map.model.descendants() + cd_map.layout.descendants()
    ):
        builder = object_to_builder.get(id(element))
        if builder is None or "id_" not in getattr(
            builder, "__dataclass_fields__", {}
        ):
            continue
        if builder.id_ != "default":
            builder.id_ = f"e{n}"
    return momapy.builder.object_from_builder(map_builder)
