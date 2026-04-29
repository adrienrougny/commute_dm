import os
import tempfile

import boltons.iterutils
import pandas
import rpy2.robjects.packages
import rpy2.robjects.vectors
import rpy2.robjects.pandas2ri

import commute_dm.queries


def make_named_gene_sets_from_collection(
    session,
    collection_name,
    namespace="ncbigene",
    with_subunits=False,
):
    query = f"""
        MATCH
            (collection:Collection)-[:HAS_ENTRY]->(entry:CollectionEntry),
            (entry)-[:HAS_ELEMENT_TO_ANNOTATIONS]->(annotations:Mapping),
            (annotations)-[:HAS_ITEM]->(item:Item),
            (item)-[:HAS_KEY]->(node)
        WHERE collection.name = '{collection_name}'
        RETURN entry AS entry, collect(node) AS nodes
    """
    result = session.execute_query(query)
    named_gene_sets = {}
    for row in result:
        entry_node = row["entry"]
        nodes = row["nodes"]
        gene_set = make_gene_set_from_nodes(
            session, nodes, namespace=namespace, with_subunits=with_subunits
        )
        named_gene_sets[entry_node["id_"]] = gene_set
    return named_gene_sets


def make_gene_set_from_nodes(
    session, nodes, namespace="ncbigene", with_subunits=False
):
    nodes = list(nodes)
    if with_subunits:
        nodes_and_subunits = commute_dm.queries.get_subunits(session, nodes)
        for _, subunits in nodes_and_subunits:
            nodes += subunits
    gene_set = set()
    identifiers = commute_dm.queries.get_identifiers(session, nodes, namespace)
    for _, identifiers in identifiers:
        for identifier in identifiers:
            gene_set.add(identifier)
    return gene_set


def make_gmt_df_from_named_gene_sets(named_gene_sets: dict):
    rows = [
        [id_, id_] + boltons.iterutils.flatten(named_gene_sets[id_])
        for id_ in named_gene_sets
    ]
    gmt_df = pandas.DataFrame(rows)
    return gmt_df


def make_gmt_file_from_named_gene_sets(named_gene_sets, output_file_path):
    gmt_df = make_gmt_df_from_named_gene_sets(named_gene_sets)
    make_gmt_file_from_gmt_df(gmt_df, output_file_path)


def make_gmt_file_from_gmt_df(gmt_df, output_file_path):
    if len(gmt_df) > 0:
        row = gmt_df.iloc[[0]]
        row = row.dropna(axis="columns")
        row.to_csv(output_file_path, sep="\t", mode="w", index=False, header=False)
        for index in range(1, len(gmt_df)):
            row = gmt_df.iloc[[index]]
            row = row.dropna(axis="columns")
            row.to_csv(output_file_path, sep="\t", mode="a", index=False, header=False)
    else:
        gmt_df.to_csv(output_file_path, sep="\t", mode="w", index=False, header=False)


def make_goat_analysis(
    gmt_df_or_file_path,
    source,
    gene_list_file_path,
    score_type="effectsize",
    p_value_cutoff=0.05,
):
    def rpy2_df_to_pandas_df(rpy2_df):
        with (
            rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter
        ).context():
            pandas_df = rpy2.robjects.conversion.get_conversion().rpy2py(rpy2_df)
        return pandas_df

    def pandas_df_to_rpy2_df(pandas_df):
        with (
            rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter
        ).context():
            rpy2_df = rpy2.robjects.conversion.get_conversion().py2rpy(pandas_df)
        return rpy2_df

    r_packages = ["goat"]
    utils = rpy2.robjects.packages.importr("utils")
    utils.chooseCRANmirror(ind=1)
    r_packages_to_install = [
        r_package
        for r_package in r_packages
        if not rpy2.robjects.packages.isinstalled(r_package)
    ]
    utils.install_packages(rpy2.robjects.vectors.StrVector(r_packages_to_install))
    goat = rpy2.robjects.packages.importr("goat")
    temp_file_path = None
    if isinstance(gmt_df_or_file_path, pandas.DataFrame):
        _, temp_file_path = tempfile.mkstemp()
        make_gmt_file_from_gmt_df(gmt_df_or_file_path, temp_file_path)
        gmt_df_or_file_path = temp_file_path
    r_gene_sets_df = goat.load_genesets_gmtfile(gmt_df_or_file_path, source)
    if temp_file_path is not None:
        os.remove(temp_file_path)
    r_gene_list_df = utils.read_csv(gene_list_file_path)
    pandas_gene_list_df = rpy2_df_to_pandas_df(r_gene_list_df)
    pandas_gene_list_df["signif"] = pandas_gene_list_df["signif"].map(
        {"True": True, "False": False}
    )
    if score_type == "effectsize_and_pvalue":
        pandas_gene_list_df["effectsize"] = pandas_gene_list_df["effectsize"] * (
            -pandas_gene_list_df["pvalue"] + 1
        )
        score_type = "effectsize"
    r_gene_list_df = pandas_df_to_rpy2_df(pandas_gene_list_df)
    r_gene_sets_filtered_df = goat.filter_genesets(
        r_gene_sets_df, r_gene_list_df, min_overlap=10, max_overlap=1500
    )
    if r_gene_sets_filtered_df.nrow > 0:
        r_result_df = goat.test_genesets(
            r_gene_sets_filtered_df,
            r_gene_list_df,
            method="goat",
            score_type=score_type,
            padj_method="BH",
            padj_cutoff=p_value_cutoff,
        )
        _, temp_file_path = tempfile.mkstemp(suffix=".csv")
        goat.save_genesets(r_result_df, r_gene_list_df, filename=temp_file_path)
        goat_df = pandas.read_csv(temp_file_path)
        os.remove(temp_file_path)
    else:
        goat_df = pandas.DataFrame()
    return goat_df
