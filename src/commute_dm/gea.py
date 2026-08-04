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
    """The `namespace` identifiers annotating the given nodes, as one set.

    Goes through `commute_dm.queries.get_annotated_nodes` rather than
    `get_subunits`, so a BEL node contributes the annotations of the entity terms
    it is built from: an `Activity` is drawn from its subject always, a `Complex`
    from its members when `with_subunits`. On the CellDesigner side that is the
    same thing `get_subunits` did -- a species plus, optionally, its subunits.
    """
    nodes = list(nodes)
    annotated_nodes = []
    for _, nodes_for_node in commute_dm.queries.get_annotated_nodes(
        session, nodes, with_subunits=with_subunits
    ):
        annotated_nodes += nodes_for_node
    gene_set = set()
    identifiers = commute_dm.queries.get_identifiers(
        session, annotated_nodes, namespace
    )
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


CRAN_REPOS = "https://cloud.r-project.org"


def _ensure_r_packages(r_packages, utils):
    """Install `r_packages` from CRAN if missing, without ever prompting.

    Two things here are what keep a notebook from hanging, and both are about
    `install.packages` rather than about the packages:

    * **It is only called when something is actually missing.**
      `install.packages` validates the library path *before* it looks at what it
      was asked to install, so calling it with an empty vector still warns
      ``'lib = "/usr/lib/R/library"' is not writable`` and then asks "Would you
      like to use a personal library instead?" -- an interactive prompt with no
      stdin behind it in a notebook, which simply blocks. It used to be called
      unconditionally, on every gene list of every level of every mode of every
      pairing.
    * **An explicit, writable `lib` is passed**, created if needed. The default
      is `.libPaths()[1]`, which is the system library on a distro R and is not
      writable by the user -- the same prompt again, this time for real. R's own
      per-user library (`R_LIBS_USER`) is the right target and is already on
      `.libPaths()`, but only once it exists.

    `repos` is passed explicitly rather than through `chooseCRANmirror`, which is
    itself interactive when it cannot pick a mirror.
    """
    missing_r_packages = [
        r_package
        for r_package in r_packages
        if not rpy2.robjects.packages.isinstalled(r_package)
    ]
    if not missing_r_packages:
        return
    base = rpy2.robjects.packages.importr("base")
    library_path = str(base.Sys_getenv("R_LIBS_USER")[0])
    if not library_path or library_path == "NA":
        raise RuntimeError(
            "R_LIBS_USER is unset, so there is no writable library to install "
            f"{missing_r_packages} into. Install them in R by hand, or set "
            "R_LIBS_USER."
        )
    base.dir_create(library_path, recursive=True, showWarnings=False)
    utils.install_packages(
        rpy2.robjects.vectors.StrVector(missing_r_packages),
        lib=library_path,
        repos=CRAN_REPOS,
    )


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

    utils = rpy2.robjects.packages.importr("utils")
    _ensure_r_packages(["goat"], utils)
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
