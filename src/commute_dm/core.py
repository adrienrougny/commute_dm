import typing
import collections
import glob
import os.path
import pathlib

import pandas
import pybiomart
import frozendict
import momapy.celldesigner
import momapy.io
import momapy.io.core
import momapy.core
import momapy.geometry
import momapy.builder
import momapy_kb.lpg.session  # noqa: F401
import momapy_kb.lpg.backends.neo4j  # noqa: F401
import commute_dm.ig
import commute_dm.queries
import commute_dm.gea  # noqa: F401

import momapy_kb.core

import commute_dm.utils  # noqa: F401


def get_interface(session):
    query = """
       CALL () {
            MATCH
                (dm_cd_collection:Collection),
                (dm_cd_collection)-[:HAS_ENTRY]->(dm_cd_entry:CollectionEntry),
                (dm_cd_entry)-[:HAS_ELEMENT_TO_ANNOTATIONS]->(dm_cd_annotations:Mapping),
                (dm_cd_entry)-[:HAS_OBJ]->(dm_cd_map:CellDesignerMap)-[:HAS_MODEL]->(dm_cd_model:CellDesignerModel),
                (dm_cd_annotations)-[:HAS_ITEM]->(dm_cd_annotations_item:Item),
                (dm_cd_annotations_item)-[:HAS_KEY]->(dm_cd_protein:Protein),
                (dm_cd_annotations_item)-[:HAS_VALUE]->(dm_cd_annotations_bag:Bag),
                (dm_cd_annotations_bag)-[:HAS_ITEM]->(dm_cd_annotation:RDFAnnotation)
            WITH
                split(dm_cd_annotation.resources[0], ":")[2] AS dm_cd_namespace,
                split(dm_cd_annotation.resources[0], ":")[-1] AS dm_cd_identifier,
                dm_cd_collection AS dm_cd_collection,
                dm_cd_entry AS dm_cd_entry,
                dm_cd_protein AS dm_cd_protein
            WHERE
                dm_cd_namespace = "hgnc.symbol"
            RETURN
                dm_cd_identifier AS identifier, dm_cd_collection AS collection, dm_cd_entry AS entry, dm_cd_protein AS protein
            UNION
            MATCH
                (ad_collection:Collection {name: "AD_KG_BEL"}),
                (ad_collection)-[:HAS_ENTRY]->(ad_entry:CollectionEntry),
                (ad_entry)-[:HAS_OBJ]->(ad_model:BELModel),
                (ad_model)-[:HAS_SUBGRAPH]->(ad_subgraph),
                (ad_subgraph)-[:HAS_NODE]->(ad_protein:Protein)
            WHERE
                ad_protein.namespace = "HGNC"
            RETURN
                ad_protein.name AS identifier, ad_collection AS collection, ad_entry AS entry, ad_protein AS protein
        }
        WITH
            identifier AS identifier,
            collect(DISTINCT [collection]) AS collections,
            collect(DISTINCT [collection, entry, protein]) AS collections_entries_proteins
        WHERE
            size(collections) >= 3
        RETURN
            identifier, collections_entries_proteins
    """
    result = session.execute_query(query)
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


def get_subgraphs_relative_to_interface(
    session,
    collection_name,
    mode: typing.Literal["downstream", "upstream"],
    max_level=-1,
    exclude_labels=None,
):
    interface = get_interface(session)
    subgraphs = collections.defaultdict(list)
    for identifier in interface:
        nodes_with_context = interface[identifier]
        for node_with_context in nodes_with_context:
            collection = node_with_context["collection"]
            if collection["name"] == collection_name:
                entry = node_with_context["entry"]
                source_node = node_with_context["node"]
                nodes, relationships = commute_dm.queries.get_subgraph(
                    session,
                    source_node,
                    relationship_types=commute_dm.ig.INFLUENCES,
                    mode=mode,
                    max_level=max_level,
                    exclude_labels=exclude_labels,
                )
                subgraphs[identifier].append(
                    {
                        "source_node_with_context": {
                            "collection": collection,
                            "entry": entry,
                            "node": source_node,
                        },
                        "nodes": nodes,
                        "relationships": relationships,
                    }
                )
    return subgraphs


def get_pd_subgraphs_downstream_of_interface(
    session, max_level=-1, exclude_labels=None
):
    return get_subgraphs_relative_to_interface(
        session,
        "PD_DM_CD",
        "downstream",
        max_level=max_level,
        exclude_labels=exclude_labels,
    )


def get_covid_subgraphs_upstream_of_interface(
    session, max_level=-1, exclude_labels=None
):
    return get_subgraphs_relative_to_interface(
        session,
        "COVID_DM_CD",
        "upstream",
        max_level=max_level,
        exclude_labels=exclude_labels,
    )


def get_n_random_nodes_from_collection(session, collection_name, n):
    query = f"""
        MATCH (collection:Collection {{name: "{collection_name}"}})-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model:Model)-[:HAS_SPECIES]->(node:Protein)
        RETURN node
        ORDER BY rand()
        LIMIT {n}
    """
    result = session.execute_query(query)
    return [row["node"] for row in result]


def make_gene_sets_from_subgraphs_relative_to_interface(session, subgraphs):
    named_gene_sets = []
    for identifier in subgraphs:
        for subgraph in subgraphs[identifier]:
            source_node_with_context = subgraph["source_node_with_context"]
            source_node = source_node_with_context["node"]
            entry = source_node_with_context["entry"]
            nodes = subgraph["nodes"]
            gene_set_id = f"{identifier}_{source_node['id_']}_{entry['id_']}"
            gene_set = commute_dm.gea.make_gene_set_from_nodes(session, nodes)
            named_gene_set = {
                "id": gene_set_id,
                "description": gene_set_id,
                "genes": gene_set,
            }
            named_gene_sets.append(named_gene_set)
    return named_gene_sets


class _SyntheticNode:
    def __init__(self, identifier):
        self.element_id = f"synthetic:{identifier}"
        self.labels = frozenset()
        self._graph = None
        self._props = {"id_": identifier, "name": identifier}

    def __getitem__(self, key):
        return self._props[key]

    def get(self, key, default=None):
        return self._props.get(key, default)


class _SyntheticRelationship:
    def __init__(self, start_node, end_node, type_):
        self.start_node = start_node
        self.end_node = end_node
        self.type = type_
        self._graph = None


def _split_interface_seeds(interface):
    seeds = {}
    for identifier, nodes_with_context in interface.items():
        covid = list(
            {
                nwc["node"]
                for nwc in nodes_with_context
                if nwc["collection"]["name"] == "COVID_DM_CD"
            }
        )
        pd = list(
            {
                nwc["node"]
                for nwc in nodes_with_context
                if nwc["collection"]["name"] == "PD_DM_CD"
            }
        )
        seeds[identifier] = {"covid": covid, "pd": pd}
    return seeds


def _make_synthetic_central_layout_element(identifier):
    layout = momapy.builder.new_builder_object(momapy.celldesigner.UnknownLayout)
    layout.position = momapy.geometry.Point(0, 0)
    layout.width = 80
    layout.height = 40
    layout.label = momapy.core.TextLayout(text=identifier, position=layout.position)
    return momapy.builder.object_from_builder(layout)


def _make_synthetic_central_model_element(identifier):
    # Pairs with the UnknownLayout above; Unknown species needs no template.
    return momapy.celldesigner.Unknown(id_=f"synthetic:{identifier}", name=identifier)


def _adjust_max_level(max_level):
    # Old (dummy) scheme: dummy at level 0, real seeds at level 1.
    # New scheme: real seeds at level 0, so subtract 1 to keep node counts equal.
    if max_level < 0:
        return max_level
    return max_level - 1


def make_and_write_cd_maps_from_interface(
    session,
    interface,
    output_dir_path,
    max_levels=None,
    min_n_nodes=None,
    covid_nodes_color="blue",
    pd_nodes_color="green",
    common_nodes_color="pink",
    interface_nodes_color="red",
):
    if max_levels is None:
        max_levels = [-1]
    commute_dm.queries.prewarm_session(session)
    seeds_by_identifier = _split_interface_seeds(interface)
    for identifier in interface:
        seeds = seeds_by_identifier[identifier]
        for max_level in max_levels:
            pd_ids = []
            covid_ids = []
            pd_and_covid_ids = []
            adjusted_max_level = _adjust_max_level(max_level)
            downstream_nodes, downstream_relationships = commute_dm.queries.get_subgraph(
                session,
                seeds["pd"],
                commute_dm.ig.INFLUENCES,
                mode="downstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            upstream_nodes, upstream_relationships = commute_dm.queries.get_subgraph(
                session,
                seeds["covid"],
                commute_dm.ig.INFLUENCES,
                mode="upstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_nodes) < min_n_nodes
                    or len(downstream_nodes) < min_n_nodes
                ):
                    continue
            nodes = list(set(downstream_nodes + upstream_nodes))
            node_to_collections = commute_dm.queries.get_collections_for_nodes(
                session, nodes
            )
            for node, collection_names in node_to_collections.items():
                node_id = node["id_"]
                if collection_names == {"COVID_DM_CD"}:
                    covid_ids.append(node_id)
                elif collection_names == {"PD_DM_CD"}:
                    pd_ids.append(node_id)
                elif collection_names == {"PD_DM_CD", "COVID_DM_CD"}:
                    pd_and_covid_ids.append(node_id)
            relationships = list(set(downstream_relationships + upstream_relationships))
            central_node = _SyntheticNode(identifier)
            synthetic_relationships = []
            for covid_seed in seeds["covid"]:
                synthetic_relationships.append(
                    _SyntheticRelationship(
                        covid_seed, central_node, commute_dm.ig.POSITIVE_INFLUENCE
                    )
                )
            for pd_seed in seeds["pd"]:
                synthetic_relationships.append(
                    _SyntheticRelationship(
                        central_node, pd_seed, commute_dm.ig.POSITIVE_INFLUENCE
                    )
                )
            ig = commute_dm.ig.make_ig_from_nodes_and_relationships(
                nodes + [central_node], relationships + synthetic_relationships
            )
            cd_map = commute_dm.ig.make_celldesigner_map_from_ig(
                session=session,
                ig=ig,
                label=f"max_level = {max_level}",
                color_node_ids=[
                    (
                        [identifier],
                        interface_nodes_color,
                    ),
                    (
                        covid_ids,
                        covid_nodes_color,
                    ),
                    (
                        pd_ids,
                        pd_nodes_color,
                    ),
                    (
                        pd_and_covid_ids,
                        common_nodes_color,
                    ),
                ],
                extra_node_layout_elements={
                    central_node: _make_synthetic_central_layout_element(identifier)
                },
                extra_node_model_elements={
                    central_node: _make_synthetic_central_model_element(identifier)
                },
            )
            output_file_path = os.path.join(
                output_dir_path, f"{identifier}_max_level_{max_level}.xml"
            )
            momapy.io.core.write(
                cd_map, output_file_path, writer="celldesigner"
            )


def make_goat_analysis_from_interface(
    session,
    interface,
    gene_lists_dir_path,
    output_dir_path,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    p_value_cutoff=0.05,
    score_type="effectsize",
):
    seeds_by_identifier = _split_interface_seeds(interface)
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
        adjusted_max_level = _adjust_max_level(max_level)
        for identifier in interface:
            seeds = seeds_by_identifier[identifier]
            upstream_nodes, _ = commute_dm.queries.get_subgraph(
                session,
                seeds["covid"],
                commute_dm.ig.INFLUENCES,
                mode="upstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            downstream_nodes, _ = commute_dm.queries.get_subgraph(
                session,
                seeds["pd"],
                commute_dm.ig.INFLUENCES,
                mode="downstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_nodes) < min_n_nodes
                    or len(downstream_nodes) < min_n_nodes
                ):
                    continue
            if mode == "downstream":
                nodes = downstream_nodes
            elif mode == "upstream":
                nodes = upstream_nodes
            elif mode == "upstream_and_downstream":
                nodes = upstream_nodes + downstream_nodes
            kept_identifiers.add(identifier)
            gene_set = commute_dm.gea.make_gene_set_from_nodes(
                session, nodes, with_subunits=with_subunits
            )
            named_gene_sets[identifier] = gene_set
        gmt_df = commute_dm.gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
        for gene_list_file_path in glob.glob(
            os.path.join(gene_lists_dir_path, "*.csv")
        ):
            gene_list_file_name = os.path.basename(gene_list_file_path)
            goat_result_df = commute_dm.gea.make_goat_analysis(
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
            for index, row in goat_result_df.iterrows():
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
    named_gene_sets = commute_dm.gea.make_named_gene_sets_from_collection(
        session, "PD_DM_CD", with_subunits=with_subunits
    )
    summary = {}
    for name in named_gene_sets:
        summary[name] = {}
    gmt_df = commute_dm.gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
    gene_list_file_names = []
    for gene_list_file_path in glob.glob(os.path.join(gene_lists_dir_path, "*.csv")):
        gene_list_file_name = os.path.basename(gene_list_file_path)
        gene_list_file_names.append(gene_list_file_name)
        goat_df = commute_dm.gea.make_goat_analysis(
            gmt_df_or_file_path=gmt_df,
            source="PD_DM_CD",
            gene_list_file_path=gene_list_file_path,
            score_type=score_type,
            p_value_cutoff=p_value_cutoff,
        )
        output_file_path = os.path.join(output_dir_path, gene_list_file_name)
        goat_df.to_csv(output_file_path)
        for index, row in goat_df.iterrows():
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
    gene_lists_dir_path,
    output_dir_path,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    min_n_hgnc=None,
):
    seeds_by_identifier = _split_interface_seeds(interface)
    if max_levels is None:
        max_levels = [-1]
    summary = collections.defaultdict(lambda: collections.defaultdict(dict))
    for max_level in max_levels:
        named_gene_sets = {}
        adjusted_max_level = _adjust_max_level(max_level)
        for identifier in interface:
            seeds = seeds_by_identifier[identifier]
            upstream_nodes, _ = commute_dm.queries.get_subgraph(
                session,
                seeds["covid"],
                commute_dm.ig.INFLUENCES,
                mode="upstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            downstream_nodes, _ = commute_dm.queries.get_subgraph(
                session,
                seeds["pd"],
                commute_dm.ig.INFLUENCES,
                mode="downstream",
                max_level=adjusted_max_level,
                filter_output_relationships=True,
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_nodes) < min_n_nodes
                    or len(downstream_nodes) < min_n_nodes
                ):
                    continue
            if mode == "downstream":
                nodes = downstream_nodes
            elif mode == "upstream":
                nodes = upstream_nodes
            elif mode == "upstream_and_downstream":
                nodes = upstream_nodes + downstream_nodes
            gene_set = commute_dm.gea.make_gene_set_from_nodes(
                session,
                nodes,
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
