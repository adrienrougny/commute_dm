import typing
import collections
import glob
import os.path
import pathlib

import pandas
import pybiomart
import frozendict
import momapy.celldesigner.core
import momapy.celldesigner.io.pickle
import momapy.io
import momapy.core
import momapy.geometry
import momapy.builder
import momapy.rendering.skia
import momapy.rendering.core
import momapy_kb.neo4j.core
import neo4j_dm.core
import neo4j_dm.ig
import neo4j_dm.queries
import neo4j_dm.gea

import commute_dm.utils


def get_interface():
    query = """
       CALL () {
            MATCH
                (dm_cd_collection:Collection),
                (dm_cd_collection)-[:HAS_ENTRY]->(dm_cd_entry:CollectionEntry),
                (dm_cd_entry)-[:HAS_RDF_ANNOTATIONS]->(dm_cd_annotations:Mapping),
                (dm_cd_entry)-[:HAS_MODEL]->(dm_cd_model:CellDesignerModel),
                (dm_cd_entry)-[:HAS_IDS]->(dm_cd_ids:Mapping),
                (dm_cd_annotations)-[:HAS_ITEM]->(dm_cd_annotations_item:Item),
                (dm_cd_annotations_item)-[:HAS_KEY]->(dm_cd_protein:Protein),
                (dm_cd_annotations_item)-[:HAS_VALUE]->(dm_cd_annotations_bag:Bag),
                (dm_cd_annotations_bag)-[:HAS_ELEMENT]->(dm_cd_annotation:RDFAnnotation)
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
                (ad_entry)-[:HAS_MODEL]->(ad_model:BELModel),
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
    result, meta = momapy_kb.neo4j.core.run(query)
    interface = {}
    for row in result:
        identifier = row[0]
        interface[identifier] = [
            {
                "collection": collection_entry_protein[0],
                "entry": collection_entry_protein[1],
                "node": collection_entry_protein[2],
            }
            for collection_entry_protein in row[1]
        ]
    return interface


def get_subgraphs_relative_to_interface(
    collection_name,
    mode: typing.Literal["downstream", "upstream"],
    max_level=-1,
    exclude_labels=None,
):
    interface = get_interface()
    subgraphs = collections.defaultdict(list)
    for identifier in interface:
        nodes_with_context = interface[identifier]
        for node_with_context in nodes_with_context:
            collection = node_with_context["collection"]
            if collection["name"] == collection_name:
                entry = node_with_context["entry"]
                source_node = node_with_context["node"]
                nodes, relationships = neo4j_dm.queries.get_subgraph(
                    source_node,
                    relationship_types=neo4j_dm.ig.INFLUENCES,
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
    max_level=-1, exclude_labels=None
):
    return get_subgraphs_relative_to_interface(
        "PD_DM_CD",
        "downstream",
        max_level=max_level,
        exclude_labels=exclude_labels,
    )


def get_covid_subgraphs_upstream_of_interface(
    max_level=-1, exclude_labels=None
):
    return get_subgraphs_relative_to_interface(
        "COVID_DM_CD",
        "upstream",
        max_level=max_level,
        exclude_labels=exclude_labels,
    )


def get_n_random_nodes_from_collection(collection_name, n):
    query = f"""
        MATCH (collection:Collection {{name: "{collection_name}"}})-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model:Model)-[:HAS_SPECIES]->(node:Protein)
        RETURN node
        ORDER BY rand()
        LIMIT {n}
    """
    result, _ = momapy_kb.neo4j.core.run(query)
    return utils.flatten_list(result)


def make_gene_sets_from_subgraphs_relative_to_interface(subgraphs):
    named_gene_sets = []
    for identifier in subgraphs:
        for subgraph in subgraphs[identifier]:
            source_node_with_context = subgraph["source_node_with_context"]
            source_node = source_node_with_context["node"]
            entry = source_node_with_context["entry"]
            nodes = subgraph["nodes"]
            gene_set_id = f"{identifier}_{source_node['id_']}_{entry['id_']}"
            gene_set = neo4j_dm.gea.make_gene_set_from_nodes(nodes)
            named_gene_set = {
                "id": gene_set_id,
                "description": gene_set_id,
                "genes": gene_set,
            }
            named_gene_sets.append(named_gene_set)
    return named_gene_sets


def _get_dummy_identifier_node(identifier):
    query = f"""
        MATCH (collection:Collection {{name: 'DUMMY_MAP'}})-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model)-[:HAS_SPECIES]->(species {{name: '{identifier}'}})
        RETURN species
    """
    result, _ = momapy_kb.neo4j.core.run(query)
    if result:
        return result[0][0]
    return None


def _make_dummy_map_from_interface(interface):
    dummy_map = momapy.celldesigner.core.CellDesignerMapBuilder()
    dummy_model = dummy_map.new_model()
    dummy_map.model = dummy_model
    dummy_layout = dummy_map.new_layout()
    dummy_map.layout = dummy_layout
    dummy_layout_model_mapping = momapy.core.LayoutModelMappingBuilder()
    dummy_map.layout_model_mapping = dummy_layout_model_mapping
    dummy_ids = {}
    for identifier in interface:
        species = dummy_model.new_element(momapy.celldesigner.core.Unknown)
        species.id_ = identifier
        species.name = identifier
        species = momapy.builder.object_from_builder(species)
        dummy_model.species.add(species)
        species_layout = dummy_layout.new_element(
            momapy.celldesigner.core.UnknownLayout
        )
        species_layout.position = momapy.geometry.Point(0, 0)
        species_layout.width = 80
        species_layout.height = 40
        species_text_layout = momapy.core.TextLayout(
            text=identifier, position=species_layout.position
        )
        species_layout.label = species_text_layout
        species_layout = momapy.builder.object_from_builder(species_layout)
        dummy_layout.layout_elements.append(species_layout)
        dummy_map.add_mapping(species, species_layout)
        dummy_ids[species] = [species.id_]
    dummy_map = momapy.builder.object_from_builder(dummy_map)
    dummy_ids = frozendict.frozendict(
        {key: frozenset(val) for key, val in dummy_ids.items()}
    )
    return dummy_map, dummy_ids


def _save_interface_to_db(interface, dummy_map, dummy_ids):
    collection_entry = neo4j_dm.core.CollectionEntry(
        id_="dummy_map",
        model=dummy_map.model,
        rdf_annotations=None,
        file_path=None,
        ids=dummy_ids,
    )
    neo4j_dm.core.save_collections_from_entries(
        [("DUMMY_MAP", [collection_entry])]
    )
    for identifier in interface:
        identifier_node = _get_dummy_identifier_node(identifier)
        pd_nodes = set(
            [
                node_with_context["node"]
                for node_with_context in interface[identifier]
                if node_with_context["collection"]["name"] == "PD_DM_CD"
            ]
        )
        for pd_node in pd_nodes:
            commute_dm.utils.merge_relationship(
                identifier_node, pd_node, neo4j_dm.ig.POSITIVE_INFLUENCE
            )
        covid_nodes = set(
            [
                node_with_context["node"]
                for node_with_context in interface[identifier]
                if node_with_context["collection"]["name"] == "COVID_DM_CD"
            ]
        )
        for covid_node in covid_nodes:
            commute_dm.utils.merge_relationship(
                covid_node, identifier_node, neo4j_dm.ig.POSITIVE_INFLUENCE
            )


def _remove_interface_from_db():
    query = """
        MATCH
            (collection:Collection {name: "DUMMY_MAP"})-[:HAS_ENTRY]->(entry)-[:HAS_MODEL]->(model)-[:HAS_SPECIES]->(species),
            (entry)-[:HAS_IDS]->(mapping)-[:HAS_ITEM]->(item)-[:HAS_VALUE]->(ids)-[:HAS_ELEMENT]->(id)
         DETACH DELETE collection, entry, model, species, mapping, item, ids, id
    """
    momapy_kb.neo4j.core.run(query)


def _make_maps_and_ids_from_dir_paths(dir_paths):
    entry_id_to_map = {}
    entry_id_to_ids = {}
    for dir_path in dir_paths:
        for input_file_path in glob.glob(os.path.join(dir_path, "*.pickle")):
            input_file_name = os.path.basename(input_file_path)
            entry_id, _ = os.path.splitext(input_file_name)
            read_result = momapy.io.read(input_file_path)
            entry_id_to_map[entry_id] = read_result.obj
            entry_id_to_ids[entry_id] = read_result.ids
    return entry_id_to_map, entry_id_to_ids


def make_and_render_igs_from_interface(
    interface,
    covid_dir_path,
    pd_dir_path,
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
    entry_id_to_map, entry_id_to_ids = _make_maps_and_ids_from_dir_paths(
        [covid_dir_path, pd_dir_path]
    )
    dummy_map, dummy_ids = _make_dummy_map_from_interface(interface)
    _save_interface_to_db(interface, dummy_map, dummy_ids)
    entry_id_to_map["dummy_map"] = dummy_map
    entry_id_to_ids["dummy_map"] = dummy_ids
    for identifier in interface:
        map_layouts = []
        identifier_node = _get_dummy_identifier_node(identifier)
        for max_level in max_levels:
            if identifier != "AGTR1" and max_level != 3:
                continue
            pd_ids = []
            covid_ids = []
            pd_and_covid_ids = []
            central_ids = []
            downstream_nodes, downstream_relationships = (
                neo4j_dm.queries.get_subgraph(
                    identifier_node,
                    neo4j_dm.ig.INFLUENCES,
                    mode="downstream",
                    max_level=max_level,
                    filter_output_relationships=True,
                )
            )
            upstream_nodes, upstream_relationships = (
                neo4j_dm.queries.get_subgraph(
                    identifier_node,
                    neo4j_dm.ig.INFLUENCES,
                    mode="upstream",
                    max_level=max_level,
                    filter_output_relationships=True,
                )
            )
            if min_n_nodes is not None:
                if (
                    len(upstream_nodes) < min_n_nodes
                    or len(downstream_nodes) < min_n_nodes
                ):
                    continue
            nodes = list(set(downstream_nodes + upstream_nodes))
            node_ids_and_context = neo4j_dm.queries.get_ids_and_context(nodes)
            for (
                node,
                ids_and_context,
            ) in node_ids_and_context:
                collection_names = set(
                    [
                        id_and_context[2]["name"]
                        for id_and_context in ids_and_context
                    ]
                )
                for id_, _, _ in ids_and_context:
                    if collection_names == set(["COVID_DM_CD"]):
                        covid_ids.append(id_)
                    elif collection_names == set(["PD_DM_CD"]):
                        pd_ids.append(id_)
                    elif collection_names == set(["PD_DM_CD", "COVID_DM_CD"]):
                        pd_and_covid_ids.append(id_)
                    elif collection_names == set(["DUMMY_MAP"]):
                        central_ids.append(id_)
            relationships = list(
                set(downstream_relationships + upstream_relationships)
            )
            ig = neo4j_dm.ig.make_ig_from_nodes_and_relationships(
                nodes, relationships
            )
            map_layout = neo4j_dm.ig.make_map_layout_from_ig(
                ig=ig,
                entry_id_to_map=entry_id_to_map,
                entry_id_to_ids=entry_id_to_ids,
                label=f"max_level = {max_level}",
                color_node_ids=[
                    (
                        central_ids,
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
            )
            map_layouts.append(map_layout)
        output_file_path = os.path.join(output_dir_path, f"{identifier}.pdf")
        if map_layouts:
            momapy.rendering.core.render_layout_elements(
                layout_elements=map_layouts,
                output_file=output_file_path,
                format_="pdf",
                renderer="skia",
                multi_pages=True,
                to_top_left=False,
            )
    _remove_interface_from_db()


def make_goat_analysis_from_interface(
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
    dummy_map, dummy_ids = _make_dummy_map_from_interface(interface)
    _save_interface_to_db(interface, dummy_map, dummy_ids)
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
            identifier_node = _get_dummy_identifier_node(identifier)
            upstream_nodes, _ = neo4j_dm.queries.get_subgraph(
                identifier_node,
                neo4j_dm.ig.INFLUENCES,
                mode="upstream",
                max_level=max_level,
                filter_output_relationships=True,
            )
            downstream_nodes, _ = neo4j_dm.queries.get_subgraph(
                identifier_node,
                neo4j_dm.ig.INFLUENCES,
                mode="downstream",
                max_level=max_level,
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
            gene_set = neo4j_dm.gea.make_gene_set_from_nodes(
                nodes, with_subunits=with_subunits
            )
            named_gene_sets[identifier] = gene_set
        gmt_df = neo4j_dm.gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
        for gene_list_file_path in glob.glob(
            os.path.join(gene_lists_dir_path, "*.csv")
        ):
            gene_list_file_name = os.path.basename(gene_list_file_path)
            goat_result_df = neo4j_dm.gea.make_goat_analysis(
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
            summary_order[identifier] += len(
                summary[identifier][gene_list_file_name]
            )
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
    _remove_interface_from_db()


def make_goat_analysis_from_pd(
    gene_lists_dir_path,
    output_dir_path,
    with_subunits=False,
    p_value_cutoff=0.05,
    score_type="effectsize",
):
    named_gene_sets = neo4j_dm.gea.make_named_gene_sets_from_collection(
        "PD_DM_CD", with_subunits=with_subunits
    )
    summary = {}
    for name in named_gene_sets:
        summary[name] = {}
    gmt_df = neo4j_dm.gea.make_gmt_df_from_named_gene_sets(named_gene_sets)
    gene_list_file_names = []
    for gene_list_file_path in glob.glob(
        os.path.join(gene_lists_dir_path, "*.csv")
    ):
        gene_list_file_name = os.path.basename(gene_list_file_path)
        gene_list_file_names.append(gene_list_file_name)
        goat_df = neo4j_dm.gea.make_goat_analysis(
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
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[
        "hsapiens_gene_ensembl"
    ]
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
    for gene_list_file_path in glob.glob(
        os.path.join(gene_lists_dir_path, "*.csv")
    ):
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
    interface,
    gene_lists_dir_path,
    output_dir_path,
    mode="downstream",
    max_levels=None,
    with_subunits=False,
    min_n_nodes=None,
    min_n_hgnc=None,
):
    dummy_map, dummy_ids = _make_dummy_map_from_interface(interface)
    _save_interface_to_db(interface, dummy_map, dummy_ids)
    if max_levels is None:
        max_levels = [-1]
    summary = collections.defaultdict(lambda: collections.defaultdict(dict))
    for max_level in max_levels:
        named_gene_sets = {}
        for identifier in interface:
            identifier_node = _get_dummy_identifier_node(identifier)
            upstream_nodes, _ = neo4j_dm.queries.get_subgraph(
                identifier_node,
                neo4j_dm.ig.INFLUENCES,
                mode="upstream",
                max_level=max_level,
                filter_output_relationships=True,
            )
            downstream_nodes, _ = neo4j_dm.queries.get_subgraph(
                identifier_node,
                neo4j_dm.ig.INFLUENCES,
                mode="downstream",
                max_level=max_level,
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
            gene_set = neo4j_dm.gea.make_gene_set_from_nodes(
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
                result_data["genes_intersection"].append(
                    list(genes_intersection)
                )
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
                summary_result = summary[identifier][gene_list_file_name][
                    max_level
                ]
                result.append(
                    f"{max_level}: {summary_result[0]}/{summary_result[1]}"
                )
                summary_order[identifier].append(
                    summary_result[0] / summary_result[1]
                )
            summary_data[gene_list_file_name].append(result)
    summary_order_df = pandas.DataFrame(
        {
            "identifier": summary_order.keys(),
            "score": [
                sum(result) / len(result) for result in summary_order.values()
            ],
        }
    )
    summary_order_df = summary_order_df.sort_values(
        by="score", ascending=False
    )
    output_summary_file_path = os.path.join(output_dir_path, "summary.csv")
    summary_df = pandas.DataFrame(summary_data)
    summary_df = summary_order_df.merge(
        summary_df, left_on="identifier", right_on="identifier", how="left"
    )
    summary_df.to_csv(output_summary_file_path)
    _remove_interface_from_db()
