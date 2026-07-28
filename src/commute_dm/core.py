import collections
import glob
import os.path
import pathlib

import pandas
import pybiomart
import momapy.celldesigner
import momapy.io
import momapy.io.core
import momapy.core
import momapy.geometry
import momapy.builder
import momapy_kb.lpg.session  # noqa: F401
import momapy_kb.lpg.backends.neo4j  # noqa: F401
import pd2af.utils
import commute_dm.ig
import commute_dm.queries

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
):
    """Split each interface entry into its upstream and downstream seed ids.

    The keys are `"upstream"` / `"downstream"` rather than the collections they
    happen to come from, so pointing the analysis at another pair of collections
    (e.g. COVID -> AD) is a parameter change, not a rewrite. Seeds are DB
    element ids, which is what the AF traversal index is keyed on.
    """
    seeds = {}
    for identifier, nodes_with_context in interface.items():
        upstream = sorted(
            {
                nwc["node"].element_id
                for nwc in nodes_with_context
                if nwc["collection"]["name"] == upstream_collection_name
            }
        )
        downstream = sorted(
            {
                nwc["node"].element_id
                for nwc in nodes_with_context
                if nwc["collection"]["name"] == downstream_collection_name
            }
        )
        seeds[identifier] = {"upstream": upstream, "downstream": downstream}
    return seeds


def _select_around_seeds(index, seeds, max_level):
    """Upstream and downstream element-id selections around one seed set.

    `max_level` means hops, plainly: the old `_adjust_max_level` compensated a
    dummy-seed scheme that no longer exists.
    """
    upstream_nodes = commute_dm.ig.select_nodes(
        index, seeds["upstream"], "upstream", max_level
    )
    downstream_nodes = commute_dm.ig.select_nodes(
        index, seeds["downstream"], "downstream", max_level
    )
    return upstream_nodes, downstream_nodes


def _make_synthetic_central_layout_element(display_name):
    layout = momapy.builder.new_builder_object(momapy.celldesigner.UnknownLayout)
    layout.position = momapy.geometry.Point(0, 0)
    layout.width = 80
    layout.height = 40
    layout.label = momapy.core.TextLayout(text=display_name, position=layout.position)
    return momapy.builder.object_from_builder(layout)


def _make_synthetic_central_model_element(identifier, display_name):
    # Pairs with the UnknownLayout above; Unknown species needs no template.
    # The name is what a reader sees (the HGNC symbol). The id_ is only unique
    # within this run -- `commute_dm.ig._renumber_ids` rewrites it before the
    # map is written, like every other id.
    return momapy.celldesigner.Unknown(
        id_=f"synthetic:{identifier}", name=display_name
    )


def _make_color_element_ids(
    index,
    species_ids,
    upstream_collection_name,
    downstream_collection_name,
    upstream_nodes_color,
    downstream_nodes_color,
    common_nodes_color,
):
    """Colour groups by **membership** in the two chosen collections.

    Membership, not set equality: with `integration_mode="hash"` a species node
    is shared by every collection that contains it (an AF collection shares its
    species with the non-AF one it was derived from), so testing
    `collection_names == {name}` matches nothing. The groups come from the
    traversal index, not `queries.get_collections_for_nodes`, whose Cypher
    traverses `(:CellDesignerMap)` and so returns `{}` for BEL nodes.
    """
    upstream_only = []
    downstream_only = []
    common = []
    for species_id in species_ids:
        collection_names = index.get_collections(species_id)
        in_upstream = upstream_collection_name in collection_names
        in_downstream = downstream_collection_name in collection_names
        if in_upstream and in_downstream:
            common.append(species_id)
        elif in_upstream:
            upstream_only.append(species_id)
        elif in_downstream:
            downstream_only.append(species_id)
    return [
        (upstream_only, upstream_nodes_color),
        (downstream_only, downstream_nodes_color),
        (common, common_nodes_color),
    ]


def make_and_write_cd_maps_from_interface(
    session,
    interface,
    index,
    output_dir_path,
    upstream_collection_name=UPSTREAM_COLLECTION_NAME,
    downstream_collection_name=DOWNSTREAM_COLLECTION_NAME,
    display_names=None,
    max_levels=None,
    min_n_nodes=None,
    include_compartment_layouts=False,
    upstream_nodes_color="blue",
    downstream_nodes_color="green",
    common_nodes_color="pink",
    interface_nodes_color="red",
):
    """Render, per interface identifier, the sub-map upstream of the identifier's
    upstream seeds and downstream of its downstream seeds.

    The sub-map is a **selection** of stored AF elements (see `commute_dm.ig`),
    joined by a synthetic central node standing for the interface identifier
    itself. `index` is the AF adjacency index, built once per run with
    `commute_dm.ig.load_af_index`.

    Output is one directory per interface protein,
    `<output_dir_path>/<display_name>/max_level_<n>.xml`.

    `include_compartment_layouts` draws the compartments: compartments are
    always in the model, but without their glyphs no box is drawn and species
    are not grouped by compartment either (`make_auto_layout` clusters on
    *mapped* compartments). Off by default, since it changes every output map.

    The interface is keyed by UniProt accession, but maps and directory names use
    the readable `display_names` (HGNC symbols, see
    :func:`get_interface_display_names`) — computed here when not supplied. The
    accession does *not* survive into the map: `_renumber_ids` rewrites every
    `id_`, the synthetic node's included. The returned stats carry both columns,
    which is the join back to the accession-keyed analyses in `4_20`.

    Returns a `DataFrame` of per-map assembly counts, one row per written map.
    Worth reading rather than discarding: `n_elements_without_glyph` counts
    selected species that are drawn *only* as complex subunits in every source
    map and so cannot stand alone in the output (they and their arcs are left
    out), and `n_extra_influences_dropped` counts seed connections to the
    central node lost the same way.
    """
    if max_levels is None:
        max_levels = [-1]
    if display_names is None:
        display_names = get_interface_display_names(session, interface)
    commute_dm.queries.prewarm_session(session)
    seeds_by_identifier = _split_interface_seeds(
        interface, upstream_collection_name, downstream_collection_name
    )
    # One hydration cache for the whole run: required for speed (momapy_kb
    # hydrates by lazy per-relationship round-trips) and for correctness (the
    # layout-model mapping is identity-keyed -- see `commute_dm.ig`).
    cache = {}
    records = []
    for identifier in interface:
        seeds = seeds_by_identifier[identifier]
        for max_level in max_levels:
            upstream_nodes, downstream_nodes = _select_around_seeds(
                index, seeds, max_level
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_nodes) < min_n_nodes
                    or len(downstream_nodes) < min_n_nodes
                ):
                    continue
            nodes = commute_dm.ig.close_over_gates(
                index, upstream_nodes | downstream_nodes
            )
            modulation_ids = commute_dm.ig.induced_modulations(index, nodes)
            species_ids = {node for node in nodes if node in index.species}
            gate_ids = {node for node in nodes if node in index.gates}
            display_name = display_names[identifier]
            central_model_element = _make_synthetic_central_model_element(
                identifier, display_name
            )
            central_layout_element = _make_synthetic_central_layout_element(
                display_name
            )
            extra_influences = [
                (
                    seed,
                    central_model_element,
                    momapy.celldesigner.PositiveInfluence,
                    momapy.celldesigner.PositiveInfluenceLayout,
                )
                for seed in seeds["upstream"]
                if seed in nodes
            ] + [
                (
                    central_model_element,
                    seed,
                    momapy.celldesigner.PositiveInfluence,
                    momapy.celldesigner.PositiveInfluenceLayout,
                )
                for seed in seeds["downstream"]
                if seed in nodes
            ]
            color_element_ids = _make_color_element_ids(
                index,
                species_ids,
                upstream_collection_name,
                downstream_collection_name,
                upstream_nodes_color,
                downstream_nodes_color,
                common_nodes_color,
            ) + [([central_model_element], interface_nodes_color)]
            cd_map, stats = commute_dm.ig.make_celldesigner_map_from_selection(
                session=session,
                cache=cache,
                element_ids=species_ids | gate_ids,
                modulation_ids=modulation_ids,
                color_element_ids=color_element_ids,
                extra_elements=[(central_model_element, central_layout_element)],
                extra_influences=extra_influences,
                with_compartment_layouts=include_compartment_layouts,
            )
            # Stored arc geometry is not reusable across a merge; `make_auto_layout`
            # repositions every glyph and rebuilds every arc's segments.
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
                    **stats,
                }
            )
    return pandas.DataFrame(records)


def _select_nodes_for_gene_set(index, seeds, max_level, mode, min_n_nodes):
    """The species nodes of one identifier's selection, or `None` if too small.

    Gates are excluded: they carry no annotation, so they cannot contribute to a
    gene set.
    """
    upstream_nodes, downstream_nodes = _select_around_seeds(index, seeds, max_level)
    if min_n_nodes is not None:
        if len(upstream_nodes) < min_n_nodes or len(downstream_nodes) < min_n_nodes:
            return None
    if mode == "downstream":
        nodes = downstream_nodes
    elif mode == "upstream":
        nodes = upstream_nodes
    elif mode == "upstream_and_downstream":
        nodes = upstream_nodes | downstream_nodes
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return {node for node in nodes if node in index.species}


def make_goat_analysis_from_interface(
    session,
    interface,
    index,
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
                index, seeds, max_level, mode, min_n_nodes
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
    index,
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
                index, seeds, max_level, mode, min_n_nodes
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
