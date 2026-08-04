import collections
import dataclasses
import glob
import os.path
import pathlib

import pandas
import pybiomart
import momapy.celldesigner
import momapy.coloring
import momapy.io
import momapy.io.core
import momapy.core
import momapy.geometry
import momapy.builder
import momapy_kb.lpg.session  # noqa: F401
import momapy_kb.lpg.backends.neo4j  # noqa: F401
import pd2af.celldesigner.building_layout
import pd2af.celldesigner.building_model
import pd2af.utils
import commute_dm.bel_submaps
import commute_dm.bel_terms
import commute_dm.queries
import commute_dm.submaps

import commute_dm.utils  # noqa: F401


# Which collections an analysis runs on is a **notebook** decision, not a library
# one: every entry point below takes the upstream and downstream collection names
# (and the interface collection names) as required arguments, and each notebook
# defines them itself. No pairing is hard-coded here.
#
# The downstream side may be an activity-flow CellDesigner collection or a BEL
# knowledge graph: `load_submap_inputs` builds the BEL side as an in-memory map (see
# `commute_dm.bel_submaps`), after which both sides are ordinary members of one
# `source_map` and nothing below distinguishes them. Which of the two a name is,
# is the one collection fact the library still needs to know:
#
# The BEL collections whose downstream side is derived from their influence-graph
# projection rather than from stored CellDesigner elements.
BEL_COLLECTION_NAMES = frozenset({"AD_KG_BEL", "PD_KG_BEL", "COVID_KG_BEL", "CBM_KG_BEL"})


def _gea():
    """Import `commute_dm.gea` lazily.

    `gea` imports `rpy2` at module scope, which fails wherever R is not usable;
    that must not make `commute_dm.core` unimportable, since only the GOAT
    entry points need it.
    """
    import commute_dm.gea

    return commute_dm.gea


def get_interface(
    session,
    collection_names,
):
    """Compute the shared-UniProt protein interface between given collections.

    All collections (CellDesigner maps and the BEL KGs alike) now carry their
    protein cross-references as uniform UniProt RDF annotations
    (`urn:miriam:uniprot:<id>` on a `:Protein` annotation key), so the interface
    is a single join over those annotations -- no CellDesigner/BEL special casing,
    no reaching inside BEL nodes, no `hgnc.symbol`. A UniProt id is in the
    interface when it is shared by **all** of `collection_names`.

    Args:
        collection_names: the collections to intersect. Passing an explicit set
            is what distinguishes, e.g., the AF collections from the non-AF ones.
            Every name must be **distinct**: the query compares
            `size(collect(DISTINCT collection.name))` against
            `size($collection_names)`, so a repeated name makes this return `{}`
            with no error.

    Returns:
        A dict mapping each shared UniProt accession to the list of
        `{"collection": <Collection node>, "entry": <CollectionEntry node>,
        "node": <Protein node>}` it was found in.
    """
    collection_names = list(collection_names)
    query = """
        MATCH
            (collection:Collection)-[:HAS_ENTRY]->(entry:CollectionEntry),
            (entry)-[:HAS_ELEMENT_TO_ANNOTATIONS]->(annotations:Mapping)-[:HAS_ITEM]->(item:Item),
            (item)-[:HAS_KEY]->(protein:Protein),
            (item)-[:HAS_VALUE]->(annotations_bag:Bag)-[:HAS_ITEM]->(annotation:RDFAnnotation)
        WHERE
            collection.name IN $collection_names
        UNWIND annotation.resources AS resource
        WITH collection, entry, protein, resource
        WHERE resource STARTS WITH "urn:miriam:uniprot:"
        WITH split(resource, ":")[-1] AS identifier, collection, entry, protein
        WITH
            identifier AS identifier,
            collect(DISTINCT collection.name) AS matched_collections,
            collect(DISTINCT [collection, entry, protein]) AS collections_entries_proteins
        WHERE
            size(matched_collections) = size($collection_names)
        RETURN
            identifier, collections_entries_proteins
    """
    result = session.execute_query(
        query,
        params={"collection_names": collection_names},
    )
    interface = {}
    for row in result:
        identifier = row["identifier"]
        interface[identifier] = [
            {
                "collection": collection_entry_protein[0],
                "entry": collection_entry_protein[1],
                "node": collection_entry_protein[2],
            }
            for collection_entry_protein in row["collections_entries_proteins"]
        ]
    return interface


def get_n_random_nodes_from_collection(session, collection_name, n):
    query = f"""
        MATCH (collection:Collection {{name: "{collection_name}"}})-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model:Model)-[:HAS_SPECIES]->(node:Protein)
        RETURN node
        ORDER BY rand()
        LIMIT {n}
    """
    result = session.execute_query(query)
    return [row["node"] for row in result]


def get_interface_display_names(session, interface, namespace="hgnc.symbol"):
    """Map each interface identifier to a readable name for the rendered maps.

    The interface is keyed by UniProt accession, which is what joins the
    collections but is unreadable in a map. Every interface protein carries an
    `hgnc.symbol` annotation, so the symbol is used instead.

    A symbol is not guaranteed unique across the interface -- `BBC3` is the
    symbol of two UniProt accessions in the current data -- so an ambiguous
    symbol is suffixed with its accession, for *both* identifiers sharing it,
    rather than letting two maps collide on one file name. An identifier with no
    symbol keeps its accession.

    When the annotations give no single symbol, a BEL node's `name` property is
    tried: it holds the bare namespace identifier (`p(HGNC:"MAPT")` -> `MAPT`).
    That fallback is no longer load-bearing -- `2_00` writes `hgnc.symbol` on the
    BEL proteins and labels them `ModelElement`, so they answer
    `commute_dm.queries.get_identifiers` like any other node, and none of the 189
    COVID x AD interface proteins falls through to its accession. It is kept for
    a node the HGNC dataset does not resolve.
    """
    identifier_to_symbol = {}
    for identifier, nodes_with_context in interface.items():
        symbols = set()
        nodes = [nwc["node"] for nwc in nodes_with_context]
        for _, identifiers in commute_dm.queries.get_identifiers(
            session, nodes, namespace
        ):
            symbols.update(identifiers)
        if len(symbols) != 1:
            symbols = {
                node["name"]
                for node in nodes
                if node.get("namespace") == "HGNC" and node.get("name")
            }
        if len(symbols) == 1:
            identifier_to_symbol[identifier] = symbols.pop()
    symbol_counts = collections.Counter(identifier_to_symbol.values())
    display_names = {}
    for identifier in interface:
        symbol = identifier_to_symbol.get(identifier)
        if symbol is None:
            display_names[identifier] = identifier
        elif symbol_counts[symbol] > 1:
            display_names[identifier] = f"{symbol}_{identifier}"
        else:
            display_names[identifier] = symbol
    return display_names


def _split_interface_seeds(
    interface,
    upstream_collection_name,
    downstream_collection_name,
    node_id_to_object=None,
    source_map=None,
    downstream_node_id_expansion=None,
):
    """Split each interface entry into its upstream and downstream seed node ids.

    The keys are `"upstream"` / `"downstream"` rather than the collections they
    happen to come from, so pointing the analysis at another pair of collections
    (e.g. COVID -> AD) is a parameter change, not a rewrite. Seeds are DB node
    ids, which is what the influence graph is keyed on.

    Passing `node_id_to_object` and `source_map` keeps only the seeds that are
    species **of the model**: an interface protein that appears solely as a
    complex subunit is not one, so it has no glyph of its own and cannot seed a
    walk. The map path filters; the gene-set path, which loads no map, does not
    -- such a seed reaches nothing there either way.

    `downstream_node_id_expansion` is a `{node_id: {node_id, ...}}` map adding, for
    each downstream seed, the other nodes standing for the same entity. It exists
    because an annotation is not always carried by the node holding the wiring: in a
    BEL KG the UniProt annotation sits on `p(HGNC:X)` while the causal edges hang off
    `act(p(HGNC:X))` (see `commute_dm.bel_submaps.load_activity_seed_expansion`). The
    expansion is applied **before** the filter, so an expanded seed still has to be a
    species of the model to survive.
    """
    filtering = node_id_to_object is not None and source_map is not None
    downstream_node_id_expansion = downstream_node_id_expansion or {}

    def keep(node_id):
        if not filtering:
            return True
        return node_id_to_object.get(node_id) in source_map.model.species

    def seed_node_ids(nodes_with_context, collection_name, expand):
        node_ids = {
            nwc["node"].element_id
            for nwc in nodes_with_context
            if nwc["collection"]["name"] == collection_name
        }
        if expand:
            node_ids = node_ids.union(
                *(downstream_node_id_expansion.get(node_id, ()) for node_id in node_ids)
            )
        return sorted(node_id for node_id in node_ids if keep(node_id))

    seeds = {}
    for identifier, nodes_with_context in interface.items():
        seeds[identifier] = {
            key: seed_node_ids(nodes_with_context, collection_name, expand)
            for key, collection_name, expand in (
                ("upstream", upstream_collection_name, False),
                ("downstream", downstream_collection_name, True),
            )
        }
    return seeds


def _select_around_seeds(influences, seeds, max_level):
    """Upstream and downstream node-id selections around one seed set.

    `max_level` means hops, plainly: the old `_adjust_max_level` compensated a
    dummy-seed scheme that no longer exists.

    One `Influences` covers both collections, as one `AfIndex` used to. The two
    AF influence graphs share no node, so this is equivalent to walking each
    collection separately -- but should a future pair of collections share
    nodes, that stops being true, and the fix is to load one `Influences` per
    collection and pass the upstream one here and the downstream one there.
    """
    return (
        influences.upstream(seeds["upstream"], max_level),
        influences.downstream(seeds["downstream"], max_level),
    )


def load_submap_inputs(
    session,
    upstream_collection_name,
    downstream_collection_name,
):
    """Everything :func:`make_and_write_submaps_from_interface` needs, once per run.

    Returns `(influences, source_map, node_id_to_object, seed_expansion,
    bel_stats)`. The
    influence graph is a few small queries; the source map is every stored map of
    the CellDesigner collections hydrated as momapy objects, which takes a couple
    of minutes and a couple of hundred megabytes. Only the map path needs it -- the
    gene-set analyses call `commute_dm.submaps.load_signed_influences` on its own.

    A collection named in `BEL_COLLECTION_NAMES` has no stored CellDesigner
    elements, so its side is built out of its influence-graph projection instead
    (~0.4 s, `commute_dm.bel_submaps`) and merged into the same `source_map`. After
    the merge a BEL node id resolves through `node_id_to_object` to a species of
    `source_map.model` exactly as a stored one does, which is why nothing below
    this function distinguishes the two.

    `seed_expansion` is empty unless the downstream side is BEL; see
    :func:`_split_interface_seeds`. `bel_stats` reports what the BEL side became
    -- species collapsed by interning, templates interned from the source map,
    complexes, subunits, badges -- and is empty when there is no BEL side.
    """
    collection_names = [upstream_collection_name, downstream_collection_name]
    bel_collection_names = [
        collection_name
        for collection_name in collection_names
        if collection_name in BEL_COLLECTION_NAMES
    ]
    cd_collection_names = [
        collection_name
        for collection_name in collection_names
        if collection_name not in BEL_COLLECTION_NAMES
    ]
    influences = commute_dm.submaps.load_signed_influences(session, cd_collection_names)
    node_id_to_object = {}
    source_map = commute_dm.submaps.load_collections_as_map(
        session, cd_collection_names, node_id_to_object
    )
    seed_expansion = {}
    bel_stats = {}
    if bel_collection_names:
        nodes, edges = commute_dm.bel_submaps.load_bel_projection(
            session, bel_collection_names
        )
        terms = commute_dm.bel_terms.load_bel_terms(session, bel_collection_names)
        # **The source map must exist first.** BEL templates are interned against
        # it (`bel_terms` invariant (a)), and the set they are interned against
        # has to be the one `make_submap_from_model_elements` re-derives -- by
        # value, from the species -- not just `model.species_templates`.
        existing_templates = frozenset(
            source_map.model.species_templates
        ) | frozenset(
            pd2af.celldesigner.building_model.collect_templates_from_species(
                source_map.model.species
            )
        )
        bel_map, species_by_node_id, bel_stats = commute_dm.bel_submaps.make_bel_map(
            nodes, terms, edges, existing_templates=existing_templates
        )
        source_map = commute_dm.bel_submaps.merge_with_source_map(source_map, bel_map)
        influences = commute_dm.bel_submaps.merge_influences(
            influences, commute_dm.bel_submaps.make_bel_influences(nodes, edges)
        )
        node_id_to_object.update(species_by_node_id)
        if downstream_collection_name in BEL_COLLECTION_NAMES:
            seed_expansion = commute_dm.bel_submaps.load_activity_seed_expansion(
                session,
                _interface_node_ids(session, collection_names, downstream_collection_name),
                projected_node_ids=nodes.keys(),
            )
    return influences, source_map, node_id_to_object, seed_expansion, bel_stats


def load_gene_set_inputs(
    session,
    upstream_collection_name,
    downstream_collection_name,
):
    """`(influences, seed_expansion)` for the gene-set analyses (a few seconds).

    What :func:`load_submap_inputs` returns minus the drawing material. The GOAT
    and intersection analyses walk node ids and read annotations off the database;
    they never hydrate a momapy object, so the couple of minutes and couple of
    hundred megabytes the `source_map` costs would buy them nothing.

    Everything else is the same, and for the same reasons: a BEL downstream side
    contributes its influence-graph projection (`commute_dm.bel_submaps`) merged
    into the same `Influences` -- the two node-id spaces are disjoint -- and the
    downstream seeds are widened to the interface proteins' activity forms,
    because BEL keeps the causal wiring on `act(p(X))` while the UniProt
    annotation that puts X in the interface sits on `p(X)`.
    """
    collection_names = [upstream_collection_name, downstream_collection_name]
    bel_collection_names = [
        collection_name
        for collection_name in collection_names
        if collection_name in BEL_COLLECTION_NAMES
    ]
    cd_collection_names = [
        collection_name
        for collection_name in collection_names
        if collection_name not in BEL_COLLECTION_NAMES
    ]
    influences = commute_dm.submaps.load_signed_influences(session, cd_collection_names)
    seed_expansion = {}
    if bel_collection_names:
        nodes, edges = commute_dm.bel_submaps.load_bel_projection(
            session, bel_collection_names
        )
        influences = commute_dm.bel_submaps.merge_influences(
            influences, commute_dm.bel_submaps.make_bel_influences(nodes, edges)
        )
        if downstream_collection_name in BEL_COLLECTION_NAMES:
            seed_expansion = commute_dm.bel_submaps.load_activity_seed_expansion(
                session,
                _interface_node_ids(session, collection_names, downstream_collection_name),
                projected_node_ids=nodes.keys(),
            )
    return influences, seed_expansion


def _interface_node_ids(session, collection_names, collection_name):
    """The nodes of one collection that carry an interface UniProt annotation.

    Only used to scope the seed expansion query to the proteins that can seed a
    walk, rather than to every protein of the knowledge graph.
    """
    interface = get_interface(session, collection_names)
    return sorted(
        {
            nwc["node"].element_id
            for nodes_with_context in interface.values()
            for nwc in nodes_with_context
            if nwc["collection"]["name"] == collection_name
        }
    )


def make_and_write_submaps_from_interface(
    session,
    interface,
    influences,
    source_map,
    node_id_to_object,
    output_dir_path,
    upstream_collection_name,
    downstream_collection_name,
    display_names=None,
    max_levels=None,
    min_n_nodes=None,
    downstream_node_id_expansion=None,
    upstream_fill=momapy.coloring.lightblue,
    downstream_fill=momapy.coloring.lightgreen,
    interface_fill=momapy.coloring.red,
):
    """Write, per interface identifier, the sub-map upstream of the identifier's
    upstream seeds and downstream of its downstream seeds.

    The sub-map holds the stored AF elements the walk selected (see
    `commute_dm.submaps`), joined by a synthetic central node standing for the
    interface identifier itself. `influences`, `source_map` and
    `node_id_to_object` come from :func:`load_submap_inputs`, once per run.

    Output is one directory per interface protein,
    `<output_dir_path>/<display_name>/max_level_<n>.xml`. Compartments are
    always drawn: without their glyphs no box appears and species are not
    grouped by compartment either, since `make_auto_layout` clusters on *mapped*
    compartments.

    A species is filled by which walk reached it: `upstream_fill` for the
    upstream selection, `downstream_fill` for the downstream one. That is the
    same partition the old colouring by collection membership produced, because
    the two AF influence graphs share no node -- and it is why the old third "in
    both collections" colour is gone: it could never fire on a *selected*
    species. (201 species nodes are in both collections, but all 201 are complex
    subunits, which no modulation points at.)

    The interface is keyed by UniProt accession, but maps and directory names use
    the readable `display_names` (HGNC symbols, see
    :func:`get_interface_display_names`) -- computed here when not supplied. The
    accession does *not* survive into the map: `renumber_ids` rewrites every
    `id_`, the synthetic node's included.

    Returns a `DataFrame` with one row per written map, holding the accession,
    the display name, the level and the element counts. The counts are what a
    change in the selection or the assembly shows up in.
    """
    if max_levels is None:
        max_levels = [-1]
    if display_names is None:
        display_names = get_interface_display_names(session, interface)
    seeds_by_identifier = _split_interface_seeds(
        interface,
        upstream_collection_name,
        downstream_collection_name,
        node_id_to_object=node_id_to_object,
        source_map=source_map,
        downstream_node_id_expansion=downstream_node_id_expansion,
    )
    records = []
    for identifier in interface:
        seeds = seeds_by_identifier[identifier]
        for max_level in max_levels:
            upstream_node_ids, downstream_node_ids = _select_around_seeds(
                influences, seeds, max_level
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_node_ids) < min_n_nodes
                    or len(downstream_node_ids) < min_n_nodes
                ):
                    continue
            node_ids = influences.close_over_gates(
                upstream_node_ids | downstream_node_ids
            )
            display_name = display_names[identifier]
            central = momapy.celldesigner.Unknown(name=display_name)
            # Fitted rather than pd2af's plain synthetic layout: an `Unknown`'s
            # default glyph is 60x30, which a display name overflows.
            central_layout_element = dataclasses.replace(
                commute_dm.bel_terms.make_fitted_synthetic_layout(central, 0),
                fill=interface_fill,
            )
            # Species only: a boolean logic gate keeps the fill it is stored with.
            fills = {
                node_id_to_object[node_id]: upstream_fill
                for node_id in influences.species_only(upstream_node_ids)
            }
            fills |= {
                node_id_to_object[node_id]: downstream_fill
                for node_id in influences.species_only(downstream_node_ids)
            }
            positive_influence = momapy.celldesigner.PositiveInfluence
            extra_influences = [
                (node_id_to_object[node_id], central, positive_influence)
                for node_id in seeds["upstream"]
                if node_id in node_ids
            ] + [
                (central, node_id_to_object[node_id], positive_influence)
                for node_id in seeds["downstream"]
                if node_id in node_ids
            ]
            cd_map = commute_dm.submaps.make_submap_from_model_elements(
                source_map,
                {node_id_to_object[node_id] for node_id in node_ids},
                fills=fills,
                extra_species=[(central, central_layout_element)],
                extra_influences=extra_influences,
            )
            # Stored arc geometry is not reusable across a merge; `make_auto_layout`
            # repositions every glyph, rebuilds every arc's segments and fits the
            # root layout.
            cd_map = pd2af.utils.make_auto_layout(cd_map)
            # One directory per interface protein, named after it; the levels of
            # one protein belong together and are what gets compared.
            protein_dir_path = os.path.join(output_dir_path, display_name)
            os.makedirs(protein_dir_path, exist_ok=True)
            output_file_path = os.path.join(
                protein_dir_path, f"max_level_{max_level}.xml"
            )
            momapy.io.core.write(cd_map, output_file_path, writer="celldesigner")
            records.append(
                {
                    "identifier": identifier,
                    "display_name": display_name,
                    "max_level": max_level,
                    # top-level only; a complex's members are counted apart
                    "n_species": len(cd_map.model.species),
                    "n_subunits": sum(
                        len(list(commute_dm.submaps.iter_subunits(species)))
                        for species in cd_map.model.species
                    ),
                    "n_modifications": sum(
                        len(getattr(species, "modifications", ()) or ())
                        for species in commute_dm.submaps.iter_species_and_subunits(
                            cd_map.model.species
                        )
                    ),
                    "n_modulations": len(cd_map.model.modulations),
                    "n_gates": len(cd_map.model.boolean_logic_gates),
                    "n_compartments": len(cd_map.model.compartments),
                    "n_templates": len(cd_map.model.species_templates),
                    "n_layout_elements": len(cd_map.layout.layout_elements),
                }
            )
    return pandas.DataFrame(records)


def _select_nodes_for_gene_set(influences, seeds, max_level, mode, min_n_nodes):
    """The species node ids of one identifier's selection, or `None` if too small.

    Gates are excluded: they carry no annotation, so they cannot contribute to a
    gene set.
    """
    upstream_node_ids, downstream_node_ids = _select_around_seeds(
        influences, seeds, max_level
    )
    if min_n_nodes is not None:
        if (
            len(upstream_node_ids) < min_n_nodes
            or len(downstream_node_ids) < min_n_nodes
        ):
            return None
    if mode == "downstream":
        node_ids = downstream_node_ids
    elif mode == "upstream":
        node_ids = upstream_node_ids
    elif mode == "upstream_and_downstream":
        node_ids = upstream_node_ids | downstream_node_ids
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return influences.species_only(node_ids)


def _add_display_name_column(df, identifier_column, display_names):
    """Insert a readable `display_name` beside a frame's identifier column.

    The analyses are keyed by UniProt accession, because that is what joins the
    collections into the interface and what is guaranteed unique -- but it is not
    what anyone reads a result table with. The HGNC symbol goes in next to it
    rather than replacing it, for the same reason
    :func:`get_interface_display_names` disambiguates a shared symbol with its
    accession: `BBC3` is the symbol of two accessions in the current data, so a
    symbol-keyed table could silently conflate two interface proteins.

    A no-op on an empty frame, which is what `make_goat_analysis` returns when
    every gene set was filtered out.
    """
    if identifier_column not in df.columns:
        return df
    display_name = df[identifier_column].map(
        lambda identifier: display_names.get(identifier, identifier)
    )
    df.insert(df.columns.get_loc(identifier_column) + 1, "display_name", display_name)
    return df


def make_goat_analysis_from_interface(
    session,
    interface,
    influences,
    gene_lists_dir_path,
    output_dir_path,
    upstream_collection_name,
    downstream_collection_name,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    downstream_node_id_expansion=None,
    display_names=None,
    p_value_cutoff=0.05,
    score_type="effectsize",
):
    """GOAT enrichment of each interface protein's selection against the gene lists.

    Gene sets are `ncbigene` identifiers, which is what `goat.test_genesets` joins
    the gene lists on (their `gene` column, written by :func:`make_goat_gene_lists`
    from the HGNC dataset's `entrez_id`).

    `downstream_node_id_expansion` comes from :func:`load_gene_set_inputs` and is
    what makes a BEL downstream side work; it is empty for a CellDesigner one.

    Every output carries a `display_name` column beside `identifier` -- see
    :func:`_add_display_name_column`.
    """
    gea = _gea()
    if display_names is None:
        display_names = get_interface_display_names(session, interface)
    seeds_by_identifier = _split_interface_seeds(
        interface,
        upstream_collection_name,
        downstream_collection_name,
        downstream_node_id_expansion=downstream_node_id_expansion,
    )
    summary = {}
    for identifier in interface:
        summary[identifier] = {}
        for (
            gene_list_file_name,
            gene_list_file_path,
        ) in commute_dm.utils.list_dir(gene_lists_dir_path):
            summary[identifier][gene_list_file_name] = []
    if max_levels is None:
        max_levels = [-1]
    kept_identifiers = set([])
    for max_level in max_levels:
        named_gene_sets = {}
        for identifier in interface:
            seeds = seeds_by_identifier[identifier]
            node_ids = _select_nodes_for_gene_set(
                influences, seeds, max_level, mode, min_n_nodes
            )
            if node_ids is None:
                continue
            kept_identifiers.add(identifier)
            gene_set = gea.make_gene_set_from_nodes(
                session,
                commute_dm.queries.get_nodes(session, node_ids),
                with_subunits=with_subunits,
            )
            named_gene_sets[identifier] = gene_set
        gmt_df = gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
        for gene_list_file_path in glob.glob(
            os.path.join(gene_lists_dir_path, "*.csv")
        ):
            gene_list_file_name = os.path.basename(gene_list_file_path)
            goat_result_df = gea.make_goat_analysis(
                gmt_df_or_file_path=gmt_df,
                source="INTERFACE",
                gene_list_file_path=gene_list_file_path,
                score_type=score_type,
                p_value_cutoff=p_value_cutoff,
            )
            output_file_name = (
                f"{pathlib.Path(gene_list_file_name).stem}_{max_level}.csv"
            )
            output_file_path = os.path.join(output_dir_path, output_file_name)
            goat_result_df = _add_display_name_column(
                goat_result_df, "id", display_names
            )
            goat_result_df.to_csv(output_file_path)
            for _, row in goat_result_df.iterrows():
                if row["signif"]:
                    summary[row["id"]][gene_list_file_name].append(max_level)
    output_summary_file_path = os.path.join(output_dir_path, "summary.csv")
    summary_data = collections.defaultdict(list)
    summary_order = collections.defaultdict(int)
    for identifier in summary:
        summary_data["identifier"].append(identifier)
        summary_data["display_name"].append(
            display_names.get(identifier, identifier)
        )
        for gene_list_file_name in summary[identifier]:
            summary_data[gene_list_file_name].append(
                summary[identifier][gene_list_file_name]
            )
            summary_order[identifier] += len(summary[identifier][gene_list_file_name])
    summary_order = [
        item[0]
        for item in sorted(
            summary_order.items(),
            key=lambda item: item[1],
            reverse=True,
        )
        if item[0] in kept_identifiers
    ]
    summary_order_df = pandas.DataFrame({"identifier": summary_order})
    summary_df = pandas.DataFrame(summary_data)
    summary_df = summary_order_df.merge(
        summary_df, left_on="identifier", right_on="identifier", how="left"
    )
    summary_df.to_csv(output_summary_file_path)


def make_goat_analysis_from_pd(
    session,
    gene_lists_dir_path,
    output_dir_path,
    with_subunits=False,
    p_value_cutoff=0.05,
    score_type="effectsize",
):
    gea = _gea()
    named_gene_sets = gea.make_named_gene_sets_from_collection(
        session, "PD_DM_CD", with_subunits=with_subunits
    )
    summary = {}
    for name in named_gene_sets:
        summary[name] = {}
    gmt_df = gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
    gene_list_file_names = []
    for gene_list_file_path in glob.glob(os.path.join(gene_lists_dir_path, "*.csv")):
        gene_list_file_name = os.path.basename(gene_list_file_path)
        gene_list_file_names.append(gene_list_file_name)
        goat_df = gea.make_goat_analysis(
            gmt_df_or_file_path=gmt_df,
            source="PD_DM_CD",
            gene_list_file_path=gene_list_file_path,
            score_type=score_type,
            p_value_cutoff=p_value_cutoff,
        )
        output_file_path = os.path.join(output_dir_path, gene_list_file_name)
        goat_df.to_csv(output_file_path)
        for _, row in goat_df.iterrows():
            if row["signif"]:
                value = row["score_type"]
            else:
                value = False
            summary[row["id"]][gene_list_file_name] = value
    summary_data = collections.defaultdict(list)
    for identifier in summary:
        summary_data["map"].append(identifier)
        for gene_list_file_name in gene_list_file_names:
            value = summary[identifier].get(gene_list_file_name, False)
            summary_data[gene_list_file_name].append(value)
    summary_df = pandas.DataFrame(summary_data)
    output_summary_file_path = os.path.join(output_dir_path, "summary.csv")
    summary_df.to_csv(output_summary_file_path)


def _get_xrefs_for_ensembl_ids(ensembl_ids, batch_size):
    server = pybiomart.Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
    start = 0
    xrefs_dfs = []
    while start < len(ensembl_ids):
        stop = start + batch_size
        xrefs_partial_df = dataset.query(
            attributes=[
                "ensembl_gene_id",
                "external_gene_name",
                "entrezgene_id",
            ],
            filters={"link_ensembl_gene_id": ensembl_ids[start:stop]},
        )
        xrefs_dfs.append(xrefs_partial_df)
        start = stop
    xrefs_df = pandas.concat(xrefs_dfs)
    return xrefs_df


# def make_goat_gene_lists(
#     gene_lists_dir_path, output_dir_path, p_value_cutoff=0.01, batch_size=400
# ):
#     for gene_list_file_path in glob.glob(
#         os.path.join(gene_lists_dir_path, "*.csv")
#     ):
#         gene_list_file_name = os.path.basename(gene_list_file_path)
#         gene_list_df = pandas.read_csv(gene_list_file_path, delimiter="\t")
#         gene_list_df = gene_list_df.rename(
#             columns={
#                 "logFC": "log2FoldChange",
#                 "P.Value": "pvalue",
#                 "adj.P.Val": "padj",
#             }
#         )
#         gene_list_df["ensembl"] = gene_list_df.index
#         ensembl_ids = list(gene_list_df["ensembl"])
#         xrefs_df = _get_xrefs_for_ensembl_ids(ensembl_ids, batch_size)
#         gene_list_df = gene_list_df.merge(
#             xrefs_df, how="outer", left_on="ensembl", right_on="Gene stable ID"
#         )
#         gene_list_df = gene_list_df.rename(
#             columns={
#                 "NCBI gene (formerly Entrezgene) ID": "gene",
#                 "Gene name": "symbol",
#                 "log2FoldChange": "effectsize",
#             }
#         )
#         gene_list_df = gene_list_df[
#             ~gene_list_df["gene"].isna()
#             & ~gene_list_df["symbol"].isna()
#             & ~gene_list_df["padj"].isna()
#         ]
#         gene_list_df["signif"] = gene_list_df["padj"] < p_value_cutoff
#         gene_list_df = gene_list_df.astype({"gene": "int64"})
#         gene_list_df = gene_list_df[
#             ~gene_list_df.duplicated(subset=["gene"])
#         ]  # TODO: take lowest pvalue?
#         output_file_path = os.path.join(output_dir_path, gene_list_file_name)
#         gene_list_df.to_csv(
#             output_file_path,
#             index=False,
#             columns=["gene", "symbol", "effectsize", "pvalue", "signif"],
#         )
#


def make_goat_gene_lists(
    gene_lists_dir_path,
    output_dir_path,
    hgnc_dataset_file_path,
    p_value_cutoff=0.01,
):
    hgnc_df = pandas.read_csv(hgnc_dataset_file_path, delimiter="\t")
    for gene_list_file_path in glob.glob(os.path.join(gene_lists_dir_path, "*.csv")):
        gene_list_file_name = os.path.basename(gene_list_file_path)
        gene_list_df = pandas.read_csv(gene_list_file_path, delimiter="\t")
        gene_list_df = gene_list_df.rename(
            columns={
                "logFC": "log2FoldChange",
                "P.Value": "pvalue",
                "adj.P.Val": "padj",
            }
        )
        gene_list_df["ensembl_gene_id"] = gene_list_df.index
        gene_list_df = gene_list_df.merge(
            hgnc_df,
            how="left",
            left_on="ensembl_gene_id",
            right_on="ensembl_gene_id",
        )
        gene_list_df = gene_list_df.rename(
            columns={
                "entrez_id": "gene",
                "log2FoldChange": "effectsize",
            }
        )
        gene_list_df = gene_list_df[
            ~gene_list_df["gene"].isna()
            & ~gene_list_df["symbol"].isna()
            & ~gene_list_df["padj"].isna()
            & ~gene_list_df["uniprot_ids"].isna()
        ]
        gene_list_df["signif"] = gene_list_df["padj"] < p_value_cutoff
        gene_list_df = gene_list_df.astype({"gene": "int64"})
        gene_list_df = gene_list_df[
            ~gene_list_df.duplicated(subset=["gene"])
        ]  # TODO: take lowest pvalue?
        output_file_path = os.path.join(output_dir_path, gene_list_file_name)
        gene_list_df.to_csv(
            output_file_path,
            index=False,
            columns=["gene", "symbol", "effectsize", "pvalue", "signif"],
        )


def make_intersection_analysis_from_interface(
    session,
    interface,
    influences,
    gene_lists_dir_path,
    output_dir_path,
    upstream_collection_name,
    downstream_collection_name,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    downstream_node_id_expansion=None,
    display_names=None,
    min_n_hgnc=None,
):
    """How much of each interface protein's selection is significantly DE.

    The naive counterpart to :func:`make_goat_analysis_from_interface`: the
    selection's gene set -- `hgnc.symbol` here, since that is what the gene lists'
    `symbol` column holds -- is intersected with the genes the gene list marks
    `signif`, rather than tested against the whole effect-size ranking. One CSV
    per gene list per level, plus a `summary.csv` ordering the interface proteins
    by the mean `|intersection| / |gene set|`.

    Every output carries a `display_name` column beside `identifier` -- see
    :func:`_add_display_name_column`.
    """
    gea = _gea()
    if display_names is None:
        display_names = get_interface_display_names(session, interface)
    seeds_by_identifier = _split_interface_seeds(
        interface,
        upstream_collection_name,
        downstream_collection_name,
        downstream_node_id_expansion=downstream_node_id_expansion,
    )
    if max_levels is None:
        max_levels = [-1]
    summary = collections.defaultdict(lambda: collections.defaultdict(dict))
    for max_level in max_levels:
        named_gene_sets = {}
        for identifier in interface:
            seeds = seeds_by_identifier[identifier]
            node_ids = _select_nodes_for_gene_set(
                influences, seeds, max_level, mode, min_n_nodes
            )
            if node_ids is None:
                continue
            gene_set = gea.make_gene_set_from_nodes(
                session,
                commute_dm.queries.get_nodes(session, node_ids),
                namespace="hgnc.symbol",
                with_subunits=with_subunits,
            )
            named_gene_sets[identifier] = gene_set
        for gene_list_file_path in glob.glob(
            os.path.join(gene_lists_dir_path, "*.csv")
        ):
            result_data = collections.defaultdict(list)
            gene_list_file_name = os.path.basename(gene_list_file_path)
            gene_list_df = pandas.read_csv(gene_list_file_path)
            genes_signif = set(gene_list_df[gene_list_df["signif"]]["symbol"])
            for identifier in named_gene_sets:
                gene_set = named_gene_sets[identifier]
                if min_n_hgnc is not None and len(gene_set) < min_n_hgnc:
                    continue
                genes_intersection = gene_set.intersection(genes_signif)
                result_data["identifier"].append(identifier)
                result_data["display_name"].append(
                    display_names.get(identifier, identifier)
                )
                result_data["gene_set"].append(list(gene_set))
                result_data["genes_signif"].append(list(genes_signif))
                result_data["genes_intersection"].append(list(genes_intersection))
                summary[identifier][gene_list_file_name][max_level] = (
                    len(genes_intersection),
                    len(gene_set),
                )
            output_file_name = (
                f"{pathlib.Path(gene_list_file_name).stem}_{max_level}.csv"
            )
            result_df = pandas.DataFrame(result_data)
            output_file_path = os.path.join(output_dir_path, output_file_name)
            result_df.to_csv(output_file_path)
    summary_data = collections.defaultdict(list)
    summary_order = collections.defaultdict(list)
    for identifier in summary:
        summary_data["identifier"].append(identifier)
        summary_data["display_name"].append(
            display_names.get(identifier, identifier)
        )
        for gene_list_file_name in summary[identifier]:
            result = []
            for max_level in summary[identifier][gene_list_file_name]:
                summary_result = summary[identifier][gene_list_file_name][max_level]
                result.append(f"{max_level}: {summary_result[0]}/{summary_result[1]}")
                summary_order[identifier].append(summary_result[0] / summary_result[1])
            summary_data[gene_list_file_name].append(result)
    summary_order_df = pandas.DataFrame(
        {
            "identifier": summary_order.keys(),
            "score": [sum(result) / len(result) for result in summary_order.values()],
        }
    )
    summary_order_df = summary_order_df.sort_values(by="score", ascending=False)
    output_summary_file_path = os.path.join(output_dir_path, "summary.csv")
    summary_df = pandas.DataFrame(summary_data)
    summary_df = summary_order_df.merge(
        summary_df, left_on="identifier", right_on="identifier", how="left"
    )
    summary_df.to_csv(output_summary_file_path)
