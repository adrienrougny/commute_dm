#!/bin/bash

INPUT_FILE_PATH=results_annotation.json
OUTPUT_FILE_PATH=results_annotation_cleaned.json
cp $INPUT_FILE_PATH $OUTPUT_FILE_PATH
sed -i -e 's/"uniprot:/"urn:miriam:uniprot:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"Uniprot:/"urn:miriam:uniprot:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"GO:/"urn:miriam:obo.go:GO%3A/g' $OUTPUT_FILE_PATH
sed -i -e 's/"go:/"urn:miriam:obo.go:GO%3A/g' $OUTPUT_FILE_PATH
sed -i -e 's/"chebi:/"urn:miriam:obo.chebi:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"CHEBI:/"urn:miriam:obo.chebi:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"mesh:/"urn:miriam:mesh:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"Mesh:/"urn:miriam:mesh:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"drugbank:/"urn:miriam:drugbank:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"interpro:/"urn:miriam:interpro:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"ncbigene:/"urn:miriam:ncbigene:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"Ncbigene:/"urn:miriam:ncbigene:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"hgnc.symbol:/"urn:miriam:hgnc.symbol:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"hgnc:/"urn:miriam:hgnc:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"pubchem.compound:/"urn:miriam:pubchem.compound:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"ncit:/"urn:miriam:ncit:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"SO:/"urn:miriam:SO:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"reactome:/"urn:miriam:reactome:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"CL:/"urn:miriam:cl:CL%3A/g' $OUTPUT_FILE_PATH
sed -i -e 's/"cl:/"urn:miriam:cl:CL%3A/g' $OUTPUT_FILE_PATH
sed -i -e 's/"BTO:/"urn:miriam:bto:BTO%3A:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"HP:/"urn:miriam:hp:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"Wikipathways:/"urn:miriam:wp:/g' $OUTPUT_FILE_PATH
sed -i -e 's/"pw:/"urn:miriam:wp:/g' $OUTPUT_FILE_PATH
sed -i -e "s/\\\u2019/'/g" $OUTPUT_FILE_PATH
