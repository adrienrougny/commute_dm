import json
import copy
import os.path
import shutil

import pandas

import commute_dm.utils


def normalize_annotations(annotations):
    new_annotations = []
    for annotation in annotations:
        if annotation:
            if not annotation.startswith("urn"):
                annotation = f"urn:miriam:{annotation}"
            annotation = annotation.replace("Ncbiprotein", "ncbiprotein")
            annotation = annotation.replace("CHEBI", "chebi")
            annotation = annotation.replace("NCBI_gene", "ncbigene")
            annotation = annotation.replace("NCBI_Gene", "ncbigene")
            annotation = annotation.replace("urn:miriam:wp", "urn:miriam:pw")
            annotation = annotation.replace("Hgnc.group", "hgnc.group")
            annotation = annotation.replace("uniprot:uniprot", "uniprot")
            annotation = annotation.replace("urn:miriam:GO:", "urn:miriam:obo.go:GO%3A")
            annotation = annotation.replace("urn:miriam:CL:", "urn:miriam:cl:CL%3A")
            annotation = annotation.replace("UNIPROT", "uniprot")
            new_annotations.append(annotation)
    return new_annotations


def get_annotation_from_annotations(collection, map, entity_id, annotations):
    for annotation in annotations:
        if (
            annotation["collection"] == collection
            and annotation["map"] == map
            and annotation["entity_id"] == entity_id
        ):
            return annotation
    return None


def get_namespaces(annotations):
    namespaces = set([])
    for annotation in annotations:
        for annotation_to_add in annotation["additional_annotations"]:
            for resource in annotation_to_add:
                namespace = ":".join(resource.split(":")[:-1])
                namespaces.add(namespace)
                if namespace == "urn:miriam":
                    print(annotation)
    return namespaces


def annotation_is_in_annotations(collection, map, entity_id, annotations):
    if (
        get_annotation_from_annotations(collection, map, entity_id, annotations)
        is not None
    ):
        return True
    return False


def remove_annotation_from_annotations(collection, map, entity_id, annotations):
    annotation = get_annotation_from_annotations(
        collection, map, entity_id, annotations
    )
    annotations.remove(annotation)


if __name__ == "__main__":
    with open("results_annotation_partly_reviewed.json") as f:
        all_annotations = json.load(f)

    missing_df = pandas.read_csv("missing_annotations.csv", sep=";")
    reviewed_annotations = []
    for i, row in missing_df.iterrows():
        row_dict = row.to_dict()
        if row_dict["Suggestion"] not in ["skip", "ok"]:
            annotations = row_dict["Suggestion"].split("+")
        else:
            if not pandas.isna(row_dict["additional_annotations"]):
                annotations = (
                    row_dict["additional_annotations"].replace(" ", "").split(",")
                )
            else:
                annotations = []
        row_dict["additional_annotations"] = [annotations]
        reviewed_annotations.append(row_dict)

    annotations_to_add = copy.deepcopy(all_annotations)
    for annotation in reviewed_annotations:
        remove_annotation_from_annotations(
            annotation["collection"],
            annotation["map"],
            annotation["entity_id"],
            annotations_to_add,
        )
        annotations_to_add.append(annotation)

    for annotation in annotations_to_add:
        new_additional_annotations = []
        for additional_annotations in annotation["additional_annotations"]:
            additional_annotations = normalize_annotations(additional_annotations)
            new_additional_annotations.append(additional_annotations)
        annotation["additional_annotations"] = new_additional_annotations

    INPUT_PD_DM_CD_DIR = "../../../build/data/pd_dm/celldesigner/"
    INPUT_COVID_DM_CD_DIR = "../../../build/data/covid_dm_annotated/celldesigner/"
    # INPUT_COVID_DM_CD_DIR = "../../../build/data/covid_dm/celldesigner/"
    OUTPUT_PD_DM_CD_DIR = "../../../build/data/pd_dm_final_annotated/celldesigner/"
    OUTPUT_COVID_DM_CD_DIR = (
        "../../../build/data/covid_dm_final_annotated/celldesigner/"
    )

    shutil.copytree(INPUT_PD_DM_CD_DIR, OUTPUT_PD_DM_CD_DIR, dirs_exist_ok=True)
    shutil.copytree(INPUT_COVID_DM_CD_DIR, OUTPUT_COVID_DM_CD_DIR, dirs_exist_ok=True)

    for annotation in annotations_to_add:
        if annotation["collection"] == "PD_DM_CD":
            dir_path = OUTPUT_PD_DM_CD_DIR
        else:
            dir_path = OUTPUT_COVID_DM_CD_DIR
        input_file_path = os.path.join(dir_path, f"{annotation['map']}.xml")
        output_file_path = input_file_path
        for resources in annotation["additional_annotations"]:
            # print(input_file_path, output_file_path, resources)
            commute_dm.utils.add_annotation_to_file(
                resources, annotation["entity_id"], input_file_path, output_file_path
            )
