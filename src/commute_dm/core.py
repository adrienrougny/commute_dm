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
import pd2af.utils
import commute_dm.queries
import commute_dm.submaps

import commute_dm.utils  # noqa: F401


# The activity-flow collections the analysis now runs on, and the interface
# they are joined with the AD BEL KG over. Referenced by string across modules
# -- change them here only.
UPSTREAM_COLLECTION_NAME = "COVID_DM_CD_AF"
DOWNSTREAM_COLLECTION_NAME = "PD_DM_CD_AF"
INTERFACE_COLLECTION_NAMES = (
    UPSTREAM_COLLECTION_NAME,
    DOWNSTREAM_COLLECTION_NAME,
    "AD_KG_BEL",
)


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
    collection_names=INTERFACE_COLLECTION_NAMES,
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
    """
    identifier_to_symbol = {}
    for identifier, nodes_with_context in interface.items():
        symbols = set()
        nodes = [nwc["node"] for nwc in nodes_with_context]
        for _, identifiers in commute_dm.queries.get_identifiers(
            session, nodes, namespace
        ):
            symbols.update(identifiers)
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
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
    node_id_to_object=None,
    source_map=None,
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
    """
    filtering = node_id_to_object is not None and source_map is not None

    def keep(node_id):
        if not filtering:
            return True
        return node_id_to_object.get(node_id) in source_map.model.species

    seeds = {}
    for identifier, nodes_with_context in interface.items():
        seeds[identifier] = {
            key: sorted(
                {
                    nwc["node"].element_id
                    for nwc in nodes_with_context
                    if nwc["collection"]["name"] == collection_name
                    and keep(nwc["node"].element_id)
                }
            )
            for key, collection_name in (
                ("upstream", upstream_collection_name),
                ("downstream", downstream_collection_name),
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
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
):
    """Everything :func:`make_and_write_submaps_from_interface` needs, once per run.

    Returns `(influences, source_map, node_id_to_object)`. The influence graph
    is two small queries; the source map is every stored map of both
    collections hydrated as momapy objects, which takes a few minutes and a
    couple of hundred megabytes. Only the map path needs it -- the gene-set
    analyses call `commute_dm.submaps.load_signed_influences` on its own.
    """
    collection_names = [upstream_collection_name, downstream_collection_name]
    influences = commute_dm.submaps.load_signed_influences(session, collection_names)
    node_id_to_object = {}
    source_map = commute_dm.submaps.load_collections_as_map(
        session, collection_names, node_id_to_object
    )
    return influences, source_map, node_id_to_object


def make_and_write_submaps_from_interface(
    session,
    interface,
    influences,
    source_map,
    node_id_to_object,
    output_dir_path,
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
    display_names=None,
    max_levels=None,
    min_n_nodes=None,
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
            central_layout_element = dataclasses.replace(
                pd2af.celldesigner.building_layout.make_synthetic_layout(central, 0),
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
                    "n_species": len(cd_map.model.species),
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


def make_goat_analysis_from_interface(
    session,
    interface,
    influences,
    gene_lists_dir_path,
    output_dir_path,
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    p_value_cutoff=0.05,
    score_type="effectsize",
):
    gea = _gea()
    seeds_by_identifier = _split_interface_seeds(
        interface, upstream_collection_name, downstream_collection_name
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
            goat_result_df.to_csv(output_file_path)
            for _, row in goat_result_df.iterrows():
                if row["signif"]:
                    summary[row["id"]][gene_list_file_name].append(max_level)
    output_summary_file_path = os.path.join(output_dir_path, "summary.csv")
    summary_data = collections.defaultdict(list)
    summary_order = collections.defaultdict(int)
    for identifier in summary:
        summary_data["identifier"].append(identifier)
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
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    min_n_hgnc=None,
):
    gea = _gea()
    seeds_by_identifier = _split_interface_seeds(
        interface, upstream_collection_name, downstream_collection_name
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
