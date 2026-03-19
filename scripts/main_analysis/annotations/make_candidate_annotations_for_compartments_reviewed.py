import os.path
import distutils.dir_util

import momapy_kb.neo4j.core
import credentials
import json

import commute_dm.utils

OUTPUT_FILE_PATH = "../../build/results/candidate_annotations_for_compartments.json"

label_to_label = {
    "ventral mesencephalic ": "ventral mesencephalon",
    "Striatal Medium Spiny Neuron D2": "striatal medium spiny neuron D2",
    "mS1 somatosensory motorized": "mS1 somatosensory motorized area",
    "callosal projection neuron-CPN": "callosal projection neuron",
    "Anterior Peduncular Area (AEP)": "anterior peduncular area",
    "Caudal Migratory System (CMS)": "caudal migratory system",
    "RADI cell (P80a)": "radiatum- and dentate-innervating cell",
    "dorsal and intermediate regions": "dorsal and intermediate hippocampus",
    "NSC/ RGL/ type 1 stem cell": "radial glia like cells",
    "MGE- neural precursor cell": "MGE-neural precursor cell",
    "MGE - neural precursor cell": "MGE-neural precursor cell",
    "basket cell (P57b)": "basket cell",
    "basket cell (P75b)": "basket cell",
    "Medial Migratory System (MMS)": "medial migratory stream",
    "corticospinal motor neuron (CSMN)": "corticospinal motor neuron",
    "DG neuroepithelium/1ry matrix": "dentate gyrus neuroepithelium",
    "corticospinal motor neuron_CSMN": "corticospinal motor neuron",
    "MGE - derived GABAergic cell": "MGE-derived GABAergic cell",
    "ventral portion": "ventral hippocampus",
    "MGE -neural precursor cell": "MGE-neural precursor cell",
    "hilus/(3ry) matrix": "hilus",
    "matrix": "mitochodrial matrix",
    "large vesicle": "astrocytic vesicle",
    "SN": "substantia nigra",
    "Peripheral Nervous System": "peripheral nervous system",
    "MP - medial pallium": "medial pallium",
    "VP - ventral pallium": "ventral pallium",
    "subpallium/ventral telencephalon (VP)": "ventral pallium",
    "Preoptic area (POa)": "preoptic area",
    "ganglionic eminence (GE) of the ventral zone (VZ)": "ganglionic eminence",
    "Lateral Ganglionic Eminence (LGE)": "lateral ganglionic eminence",
    "dorsal telencephalon / pallium": "dorsal telencephalon",
    "Caudal Ganglionic Eminence (CGE)": "caudal ganglionic eminence",
    "subgranular zone (SGZ)": "subgranular zone",
    "Basomedial Amygdala Nuclei (BMA)": "basomedial amygdala nuclei",
    "DP - Dorsal Pallium": "dorsal pallium",
    "LMS-lateral migratory system": "lateral migratory system",
    "RPC -rostral patterning center": "rostral patterning center",
    "VZ - ventricular zone": "venticular zone",
}

label_to_iris = {
    "cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000000"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "": {
        "source": "manual",
        "iris": [[]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "axon": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0030424"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "pancreas": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0000988"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "neural precursor cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A4042021"]],
        "set_operation": None,
        "exact": False,
    },
    "interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000099"]],
        "set_operation": None,
    },
    "cell membrane": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005886"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "extracellular space": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005576"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "brain": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0000142"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cortex": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D002540"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "ribosome": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005840"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "postsynaptic terminal": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0098794"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "human host": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:taxonomy:9606",
                "http://purl.obolibrary.org/obo/OHMI_0000485",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find exact term so put 'human' and 'host'",
    },
    "blood": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:UBERON:0000178",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "glutamatergic pyramidal neuron": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000679",
                "urn:miriam:obo.clo:CL%3A0000598",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find exact term so put 'glutamatergic neuron' and 'pyramidal neuron'",
    },
    "epithelial cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000066",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "mitochondrial nucleoid": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0042645",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "apoptotic neuron": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000540",
                "urn:miriam:obo.go:GO%3A0051402",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put 'neuron', 'particiapates in', 'neuron apoptosis process'",
    },
    "activated microglia m2": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0002629",
                "urn:miriam:obo.clo:CL%3A0000890",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "To check. Can't find the exact term, so put an intersection 'mature microglial cell' and 'alternatively activated macrophage'",
    },
    "cortical gabaergic  interneuron": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0010011",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": "Typo: two spaces between GABAergic and interneuron",
    },
    "neocortical interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:CL:0008031", "urn:miriam:mesh:D019579"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find exact term (there's only neocortical basket cell which is more specialized?), so put intersection of 'interneuron', 'located in', and 'neocortex'",
    },
    "layer ii": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C33137"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "neurogliaform cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000693"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "gabaergic interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0011005"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "layer iv": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C32844"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "pyramidal cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000598"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cortical gabaergic interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0010011"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "spinal cord": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0001279"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "stratum radiatum": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D056547"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "layer v": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C32653"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "corticothalamic neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A4023013"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put 'neuron', 'has quality', and 'corticothalamic projecting'",
    },
    "somatosensory cortex": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0004353"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "dorsal hippocampal ca1 area": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45874", "urn:miriam:mesh:D056547"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "preoptic area (poa)": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D011301"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "caudal ganglionic eminence (cge)": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C73851", "urn:miriam:mesh:D000097803"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "lateral ganglionic eminence (lge)": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C25230", "urn:miriam:mesh:D000097803"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "stratum pyramidale": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D056547", "urn:miriam:obo.bto:BTO%3A0001066"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "peripheral nervous system": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0001028"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "ivy cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A4042013"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "corticospinal motor neuron_csmn": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000100",
                "urn:miriam:mesh:D011712",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find the exact term, so put 'motor neuron', 'located in', and 'corticospinal tract'",
    },
    "retrosplenial cortex": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D006179"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cerebellum": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D002531"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "subgranular zone (sgz)": {
        "source": "manual",
        "iris": [["urn:miriam:wikipedia.en:Subgranular_zone"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "basket cell (p57b)": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000118"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check, don't know what the p57b refers to",
    },
    "basket cell (p75b)": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000118"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check, don't know what the p75b refers to. We also have a label with 'P57b' elsewhere, is one of the two a typo?",
    },
    "medial ganglionic eminence (mge)": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C25232", "urn:miriam:mesh:D000097803"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "pyramidal neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000598"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "liver": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D008099"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "motor cortex": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D009044"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "granule cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000120"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "layer iii": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C32571"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "amygdala": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D000679"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "layer i": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C33137"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "layer vi": {
        "source": "manual",
        "iris": [["urn:miriam:UBERON:0005395"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "hippocampus": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D006624"]],
        "set_operation": None,
        "exact": True,
        "comment": "Interestingly (?) there is no UBERON term for hippocamus, only for 'hippocampal formation', which is broader",
    },
    "hypothalamus": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D007031"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "schwann cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0002573"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "radi cell (p80a)": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D056547",
                "urn:miriam:mesh:D018891",
                "urn:miriam:obo.clo:CL%3A0000099",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find",
    },
    "ca1": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D056547"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "ca1 area": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D056547"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "corpus callosum": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0000615"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "subiculum": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C33648"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "glutaminergic neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000679"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check: is it indeed 'glutaminergic neuron' or should it rather be 'glutamatergic neuron'? There is literature on the glutamate-glutamine cycle but not much on glutaminergic neurons. If it is indeed 'glutaminergic neuron' that is meant, I could not find an exact match. Closest could be 'mitral cell' (CL:1001502)",
    },
    "cholinergic interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000108", "urn:miriam:CL:0000099"]],
        "set_operation": "intersection",
        "exact": True,
        "comment": "To check. Could not find an exact term so put 'cholinergic neuron' and 'interneuron'",
    },
    "nucelus": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005634"]],
        "set_operation": None,
        "exact": True,
        "comment": "Typo: nucelus->nucleus",
    },
    "gabaergic neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000617"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cholinergic neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000108"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "striatal neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0004225"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "substantia nigra": {
        "source": "manual",
        "iris": [["urn:miriam:UBERON:0002038"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "purkinje cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000121"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "vesicle": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0031982"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "post synapse": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0098794"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "large vesicle": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0070382",
                "urn:miriam:obo.clo:CL%3A0002604",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "Could not find the exact term, put 'dense core vesicle' with exact = False, not sure if it corresponds in this context. The GO has a 'dense core granule' (obo.go:GO%3A0031045) term that could be a synonym. If it is the case it is preferable to use the GO term?",
    },
    "striatum": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0001311"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "substantia nigra and striatum": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.bto:BTO%3A0000143",
                "urn:miriam:obo.bto:BTO%3A0001311",
            ]
        ],
        "set_operation": "union",
        "exact": False,
        "comment": None,
    },
    "da neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000700"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "pulmonary endothelial cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A1001567"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cholinergic cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000108"]],
        "set_operation": None,
        "exact": False,
        "comment": "To check. Could not find the exact term, so put 'cholinergic neuron'. Not sure if there are cholinergic cells that are not neurons in practice.",
    },
    "pallial-derived cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C45874",
                "urn:miriam:mesh:D013687",
                "urn:miriam:obo.clo:CL%3A0000000",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'cell', 'derived from anatomical part', and 'pallium'",
    },
    "activated microglia": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0002629"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "vp - ventral pallium": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45875", "urn:miriam:mesh:D013687"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "vz - ventricular zone": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0003654"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "anterior peduncular area (aep)": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25232",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:ncit:C45875",
                "urn:miriam:ncit:C73851",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "To check. I could not find the exact term. I put 'ansa peduncularis' because it is the only possibly suitable term which contains 'peduncular' I found (and in my understanding it is an anterior region of the thalamus)",
    },
    "dorsal cortex": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45874", "urn:miriam:mesh:D002540"]],
        "set_operation": "intersection",
        "exact": True,
        "comment": "Could not find the exact term, so put intersection of 'cerebral cortex' and 'dorsal part'",
    },
    "somatosensory cortex s1": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0004353"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cajal‐retzius cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000695"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "dentate gyrus (dg)": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0002496"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "medial amygdala": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D066276"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "basomedial amygdala nuclei (bma)": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D066276"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "mp - medial pallium": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D013687", "urn:miriam:ncit:C25232"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "medial migratory system (mms)": {
        "source": "manual",
        "iris": [["urn:miriam:pubmed:27034423"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find a suitable term. Can found terms for 'migratory stream', but no medial one.",
    },
    "mitochondrium": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005739"]],
        "set_operation": None,
        "exact": True,
        "comment": "Typo: The label is in German? 'mitochondrium'->'mitochondrion' or 'mitochondria'",
    },
    "pons": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0001101"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "mitochondrial space matrix": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005759"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "actin cytoskeleton": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0015629"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "stratum lacunosum-moleculare": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D056654"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "globus pallidus pars interna": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0002248"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "presynaptic terminal": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D017729"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "penduncolo pontine tegmental nucleus (pptg)": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D045042"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check",
    },
    "mitochondria-associated endoplasmic reticulum membrane": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0044233"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cath cell line": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0005820"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check, found term 'CATH.a cell', which is an immortal mouse cell line",
    },
    "mge- neural precursor cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25232",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:obo.clo:CL%3A4042021",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'precursor cell', 'neural cell', 'located in', and 'medial ganglionic eminence'",
    },
    "mge - neural precursor cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25232",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:obo.clo:CL%3A4042021",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'precursor cell', 'neural cell', 'located in', and 'medial ganglionic eminence'",
    },
    "substantia nigra pars reticulata": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0003750"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "glutamatergic projection": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0043005",
                "urn:miriam:obo.clo:CL%3A0000679",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find the exact term, so put intersection of 'neuron projection', 'part of', and 'glutamatergic neuron'",
    },
    "globus pallidus pars externa": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0002247"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "callosal projection neuron-cpn": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3ACL:0000540",
                "urn:miriam:obo.go:GO%3A0043005",
                "urn:miriam:mesh:D003337",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find the exact term, so put intersection of 'neuron', 'has part', 'neuron projection', 'located in', 'corpus callosum'",
    },
    "ventral mesencephalic ": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D013681"]],
        "set_operation": None,
        "exact": True,
        "comment": "Typo: remove the ending space",
    },
    "granule neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000120"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "matrix": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005759"]],
        "set_operation": None,
        "exact": True,
        "comment": "Could not find a generic 'matrix' term. Need to check the context to see what it is a matrix of",
    },
    "trans-golgi network": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0005802"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "axonal varicosity": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0043196"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "arteriole": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D001160"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "cortical neuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0010012"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "caudal migratory system (cms)": {
        "source": "manual",
        "iris": [["urn:miriam:pubmed:16079409"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find a suitable term. Can found terms for 'migratory stream', but only a specialized muscle one for caudal (posterior) migratory stream",
    },
    "ventromedial hypothalamus": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D014697"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "corticospinal motor neuron (csmn)": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000100",
                "urn:miriam:mesh:D011712",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "To check. Could not find the exact term, so put 'motor neuron', 'located in', and 'corticospinal tract'",
    },
    "ca fields": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:BTO:0003705"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "phagocytic vesicle": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0045335"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "developing cerebral cortex": {
        "source": "manual",
        "iris": [["urn:miriam:obo.go:GO%3A0007420", "urn:miriam:mesh:D002540"]],
        "set_operation": None,
        "exact": False,
        "comment": None,
    },
    "dp - dorsal pallium": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45874", "urn:miriam:mesh:D013687"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "subpallium/ventral telencephalon (vp)": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45875", "urn:miriam:mesh:D013687"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "hapi microglia": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0003618"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "sn": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0000143"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check within the context to see if it indeed stands for 'substantia nigra'",
    },
    "neuron/gt1-7": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0002708"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "rostrocaudal axis": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C94393",
                "urn:miriam:ncit:C73851",
                "urn:miriam:ncit:C25154",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find the exact term, put 'anterior-posterior axis', which is tagged as being broader than 'rostrocaudal axis'",
    },
    "dorsomedial pallium": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C45874",
                "urn:miriam:mesh:D013687",
                "urn:miriam:ncit:C25232",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'dorsal pallium' and 'medial pallium'",
    },
    "svz": {
        "source": "manual",
        "iris": [["urn:miriam:obo.bto:BTO%3A0003090"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "posterior subpallial domain": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C45875",
                "urn:miriam:mesh:D013687",
                "urn:miriam:ncit:C25622",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find the term so put intersection of 'subpallium' and 'dorsal part'",
    },
    "cortical interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0008031"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "hippocampal interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A1001569"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "dorsal telencephalon / pallium": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45874", "urn:miriam:mesh:D013687"]],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "mge - derived gabaergic cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000617",
                "urn:miriam:ncit:C25232",
                "urn:miriam:mesh:D000097803",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'gabaergic neuron', 'derived from anatomical part', and 'medial ganglionic eminence'",
    },
    "poa-derived cell": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D011301", "urn:miriam:obo.clo:CL%3A0000000"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'cell', 'derived from anatomical part', and 'preoptic area'",
    },
    "poa-derived gabaergic interneuron": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D011301", "urn:miriam:obo.clo:CL%3A0011005"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'GABAergic interneuron', 'derived from anatomical part', and 'preoptic area'",
    },
    "cge-derived gabaergic cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000617",
                "urn:miriam:ncit:C73851",
                "urn:miriam:mesh:D000097803",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'GABAergic neuron', 'derived from anatomical part', and 'caudal ganglionic eminence'",
    },
    "cge-derived cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000000",
                "urn:miriam:ncit:C73851",
                "urn:miriam:mesh:D000097803",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'cell', 'derived from anatomical part', and 'caudal ganglionic eminence'",
    },
    "lge-derived gabaergic cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25230",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:obo.clo:CL%3A0000617",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'GABAergic neuron', 'derived from anatomical part', and 'lateral ganglionic eminence'",
    },
    "poa-derived gabaergic interneuron gad1+": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D011301",
                "urn:miriam:obo.clo:CL%3A0011005",
                "urn:miriam:uniprot:Q99259",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'GABAergic interneuron', 'derived from anatomical part', 'preoptic area', 'has quality', 'positive' and 'GAD1'",
    },
    "mge -neural precursor cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25232",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:obo.clo:CL%3A4042021",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'precursor cell', 'neural cell', 'neural stem cell', 'located in', and 'medial ganglionic eminence'",
    },
    "lge-neural precursor cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C25230",
                "urn:miriam:mesh:D000097803",
                "urn:miriam:obo.clo:CL%3A4042021",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put intersection of 'precursor cell', 'neural cell', 'neural stem cell', 'located in', and 'lateral ganglionic eminence'",
    },
    "dorsal entorhinal cortex": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45874", "urn:miriam:mesh:D018728"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term so put intersection of 'entorhinal cortex' and 'dorsal part'",
    },
    "dorsal poa (poa1)": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D011301", "urn:miriam:ncit:C45874"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term so put intersection of 'entorhinal cortex' and 'dorsal part'",
    },
    "caudal cortex": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C73851", "urn:miriam:ncit:C12443"]],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find the exact term so put intersection of 'inferior side' and 'cerebral cortex'",
    },
    "ganglionic eminence (ge) of the ventral zone (vz)": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D000097803",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "s1 parietal / somatosensory cortex": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:UBERON:0008930",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "ms1 somatosensory motorized": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D013003"]],
        "set_operation": None,
        "exact": True,
        "comment": "Could not figure exactly what it is, a part of the 'somatosensory cortex'?",
    },
    "neurocortical interneuron": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0008031",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "dg neuroepithelium/1ry matrix": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D054258urn:miriam:obo.go:GO%3A0021542",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find the exact term, so put 'developing neuroepithelium', 'participates in', 'dentate gyrus development'",
    },
    "lms-lateral migratory system": {
        "source": "manual",
        "iris": [["urn:miriam:pubmed:17093077"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find a suitable term. Can found terms for 'migratory stream', but no lateral one.",
    },
    "p33b neurogliaform cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000693",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find exact term, so put broader term 'neurogliaform cell'",
    },
    "ventro-caudal poa (poa2)": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D011301",
                "urn:miriam:ncit:C73851",
                "urn:miriam:ncit:C45875",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find the term so put intersection of 'ventral side', 'posterior side' and 'preoptic area'",
    },
    "hilus/(3ry) matrix": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D018891", "urn:miriam:obo.go:GO%3A0021542"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find a corresponding term. Seems to be part of the dentate gyrus, but is it the 'hilus of dentate gyrus'? Does not seem so, probably a part only during the development?",
    },
    "nsc/ rgl/ type 1 stem cell": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A0000681"]],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find any term for radial glia-like cell or type 1 stem cell, so put 'neural stem cell', 'located in', and 'dentate gyrus subgranular zone'",
    },
    "interneuron-selective cell": {
        "source": "manual",
        "iris": [["urn:miriam:pubmed:20130170"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find a suitable term. Could be one under 'interneuron' but I'm not sure.",
    },
    "developing basal ganglia": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D001479", "urn:miriam:obo.go:GO%3A0061548"]],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find the exact term, so put 'collection of basal ganglia', 'participates in', 'ganglion development'",
    },
    "radiatum-retrohippocampal cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:mesh:D056547",
                "urn:miriam:obo.go:GO%3A0043005",
                "urn:miriam:mesh:D006179",
            ]
        ],
        "set_operation": None,
        "exact": False,
        "comment": "Can't find the exact term, put 'CNS long range interneuron' which is broader",
    },
    "rpc -rostral patterning center": {
        "source": "manual",
        "iris": [["urn:miriam:pubmed:25889070"]],
        "set_operation": None,
        "exact": True,
        "comment": "Can't find the term and not sure what it refers to",
    },
    "p5": {
        "source": "manual",
        "iris": [[]],
        "set_operation": None,
        "exact": True,
        "comment": "Protein disulfide isomerase P5?",
    },
    "ventral portion": {
        "source": "manual",
        "iris": [["urn:miriam:ncit:C45875", "urn:miriam:mesh:D006624"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check, probably looking at the context is needed",
    },
    "striatal medium spiny neuron d2": {
        "source": "manual",
        "iris": [["urn:miriam:obo.clo:CL%3A4023029"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check, don't know what the d2 refers to",
    },
    "oxidized lipid membrane": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0016020",
                "urn:miriam:sio:SIO_000292",
                "http://purl.obolibrary.org/obo/MOP_0000568",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the term, so put 'membrane', 'is target in', 'oxidation'",
    },
    "culture of neurons and glia": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.bto:BTO%3A0000214",
                "urn:miriam:obo.clo:CL%3A0000125",
                "urn:miriam:obo.clo:CL%3A0000540",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Could not find the exact term. Put 'cell culture', 'contains', 'glical cell' and 'neuron'",
    },
    "innate immune cell": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.clo:CL%3A0000000",
                "urn:miriam:obo.go:GO%3A0045087",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put 'cell', 'participates in', 'innate immune response'. An alternative would be to have an union of all innate immune cell types",
    },
    "nigrostriatal da projection": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0043005",
                "urn:miriam:obo.bto:BTO%3A0004776",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put 'neuron projection', 'part of', and 'striatonigral neuron'",
    },
    "microglia enriched culture": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.bto:BTO%3A0000214",
                "urn:miriam:obo.clo:CL%3A0000129",
            ]
        ],
        "set_operation": "intersection",
        "exact": False,
        "comment": "Can't find the exact term, so put 'cell culture', 'contains', and 'microglial cell'",
    },
    "multiform layer": {
        "source": "manual",
        "iris": [["urn:miriam:mesh:D019579"]],
        "set_operation": None,
        "exact": True,
        "comment": "To check.",
    },
    "mitochondrial derived vesicle": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:obo.go:GO%3A0099073",
            ]
        ],
        "set_operation": None,
        "exact": True,
        "comment": None,
    },
    "dorsal and intermediate regions": {
        "source": "manual",
        "iris": [
            [
                "urn:miriam:ncit:C45874",
                "urn:miriam:ncit:C73705",
                "urn:miriam:mesh:D006624",
            ]
        ],
        "set_operation": "union",
        "exact": False,
        "comment": "Can't find the exact term, so put 'dorsal part (qualifier)' and 'medial part (qualifier)'. Could perhaps be improved by checking the context",
    },
}

if __name__ == "__main__":
    # query = """
    #     MATCH
    #         (collection:Collection)-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model:CellDesignerModel),
    #         (model)-[:HAS_COMPARTMENT]->(compartment:Compartment),
    #         (entry)-[:HAS_IDS]->(ids:Mapping),
    #         (entry)-[:HAS_RDF_ANNOTATIONS]->(annotations:Mapping),
    #         (ids)-[:HAS_ITEM]->(model_id_item:Item)-[:HAS_KEY]->(model),
    #         (model_id_item)-[:HAS_VALUE]->(model_id_bag:Bag)-[:HAS_ELEMENT]->(model_id:String),
    #         (ids)-[:HAS_ITEM]->(compartment_id_item:Item)-[:HAS_KEY]->(compartment),
    #         (compartment_id_item)-[:HAS_VALUE]->(compartment_id_bag:Bag)-[:HAS_ELEMENT]->(compartment_id:String)
    #     OPTIONAL MATCH
    #         (annotations)-[:HAS_ITEM]->(compartment_annotation_item:Item)-[:HAS_KEY]->(compartment),
    #         (compartment_annotation_item)-[:HAS_VALUE]->(compartment_annotation_bag:Bag)-[:HAS_ELEMENT]->(compartment_annotation)
    #     UNWIND
    #         CASE
    #             WHEN compartment_annotation.resources = [] THEN [NULL]
    #             WHEN compartment_annotation IS NULL THEN [NULL]
    #             ELSE compartment_annotation.resources
    #         END
    #         AS compartment_annotation_resource
    #     WITH
    #         collection AS collection,
    #         model_id AS model_id,
    #         compartment_id AS compartment_id,
    #         compartment AS compartment,
    #         COLLECT(compartment_annotation_resource) AS compartment_annotation_resources
    #     RETURN
    #         collection.name, model_id.value, compartment_id.value, compartment.name,
    #         [compartment_annotation_resource IN compartment_annotation_resources
    #             WHERE NOT compartment_annotation_resource CONTAINS "pubmed"
    #             AND NOT compartment_annotation_resource CONTAINS "doi"
    #         ] AS compartment_annotation_resources
    # """
    #
    query = """
        MATCH
            (collection:Collection)-[:HAS_ENTRY]->(entry:CollectionEntry)-[:HAS_MODEL]->(model:CellDesignerModel),
            (model)-[:HAS_COMPARTMENT]->(compartment:Compartment),
            (entry)-[:HAS_IDS]->(ids:Mapping),
            (entry)-[:HAS_RDF_ANNOTATIONS]->(annotations:Mapping),
            (ids)-[:HAS_ITEM]->(compartment_id_item:Item)-[:HAS_KEY]->(compartment),
            (compartment_id_item)-[:HAS_VALUE]->(compartment_id_bag:Bag)-[:HAS_ELEMENT]->(compartment_id:String)
        OPTIONAL MATCH
            (annotations)-[:HAS_ITEM]->(compartment_annotation_item:Item)-[:HAS_KEY]->(compartment),
            (compartment_annotation_item)-[:HAS_VALUE]->(compartment_annotation_bag:Bag)-[:HAS_ELEMENT]->(compartment_annotation)
        UNWIND
            CASE
                WHEN compartment_annotation.resources = [] THEN [NULL]
                WHEN compartment_annotation IS NULL THEN [NULL]
                ELSE compartment_annotation.resources
            END
            AS compartment_annotation_resource
        WITH
            collection AS collection,
            entry.id_ AS model_id,
            compartment_id AS compartment_id,
            compartment AS compartment,
            COLLECT(compartment_annotation_resource) AS compartment_annotation_resources
        RETURN
            collection.name, model_id, compartment_id.value, compartment.name,
            [compartment_annotation_resource IN compartment_annotation_resources
                WHERE NOT compartment_annotation_resource CONTAINS "pubmed"
                AND NOT compartment_annotation_resource CONTAINS "doi"
            ] AS compartment_annotation_resources
        """
    momapy_kb.neo4j.core.connect(
        credentials.NEO4J_URI,
        credentials.NEO4J_USERNAME,
        credentials.NEO4J_PASSWORD,
    )
    collection_names_model_ids_entity_ids_entity_labels = []
    results, meta = momapy_kb.neo4j.core.run(query)
    labels_no_iris = []
    n = 0
    for result in results:
        label = result[3]
        if label is not None:
            label_lower = label.lower()
            print(label_lower)
            iris = result[4]
            if iris:
                if label_lower in label_to_iris:
                    if iris not in label_to_iris[label_lower]["iris"]:
                        label_to_iris[label_lower]["iris"].append(iris)
                else:
                    label_to_iris[label_lower] = {
                        "source": "from_db",
                        "iris": [iris],
                        "set_operation": None,
                        "exact": False,
                        "comment": None,
                    }
            else:
                n += 1
                labels_no_iris.append(label_lower)
                collection_names_model_ids_entity_ids_entity_labels.append(
                    {
                        "collection_name": result[0],
                        "model_id": result[1],
                        "entity_id": result[2],
                        "entity_label": result[3],
                    }
                )
    for label in labels_no_iris:
        if label not in label_to_iris:
            print(label)
    annotations = []
    for (
        collection_name_model_id_entity_id_entity_label
    ) in collection_names_model_ids_entity_ids_entity_labels:
        label = collection_name_model_id_entity_id_entity_label["entity_label"].lower()
        if label in label_to_iris and label:
            source_iris_set_operation_exact_comment = label_to_iris[label]
            annotation = (
                collection_name_model_id_entity_id_entity_label
                | source_iris_set_operation_exact_comment
            )
            annotations.append(annotation)
    # with open(OUTPUT_FILE_PATH, "w") as f:
    #     description = "This file was automatically generated by the script located at '/src/commute_dm_develop/make_candidate_annotations_for_compartments.py'. It gives candidate annotations for compartments in the PD and COVID disease maps that are lacking annotations. For such a compartment, its candidate annotations were either generated automatically by using existing annotations of compartments sharing the same lower cased label (indicated by 'source': 'from_db'), or were manually curated ('source': 'manual') when it didn't share any label. The annotations are given by the 'iris' list. Annotations from the list are alternatives. An annotation may contain several IRIs. In this case, its the IRIs taken altogether that constitute the annotation, and the 'set_operation' attribute indicate whether the different IRIs have to be taken as a sort of 'union' (e.g., the compartment belongs to class A or to class B) or as a sort of 'intersection' (e.g., the compartment belongs to class A and to class B). The 'exact' attribute is set to 'true' for manually curated annotations when the annotation is formed of a unique IRI that corresponds exactly to the represented concept. It is set to 'false' otherwise, in particular for annotations generated from the DB, that were not checked manually. Some annotations should be checked. They should all contain 'To check' in their 'comment'. There are also a few typos in labels, indicated by 'Typo' in the 'comment'."
    #     json.dump(
    #         {
    #             "description": description,
    #             "annotations": annotations,
    #         },
    #         f,
    #     )
    # PD_DM_INPUT_DIR = "../../data/pd_dm_sources/celldesigner/"
    COVID_DM_INPUT_DIR = "../../../build/data/covid_dm/celldesigner/"
    # PD_DM_OUTPUT_DIR = "../../build/data/pd_dm_sources_annotated/celldesigner/"
    COVID_DM_OUTPUT_DIR = "../../../build/data/covid_dm_annotated/celldesigner/"
    # distutils.dir_util.copy_tree(PD_DM_INPUT_DIR, PD_DM_OUTPUT_DIR)
    distutils.dir_util.copy_tree(COVID_DM_INPUT_DIR, COVID_DM_OUTPUT_DIR)
    for annotation in annotations:
        if annotation["collection_name"] == "PD_DM_CD":
            continue
        else:
            input_dir = COVID_DM_OUTPUT_DIR
        output_dir = input_dir
        file_name = f"{annotation['model_id']}.xml"
        input_file_path = os.path.join(input_dir, file_name)
        output_file_path = os.path.join(output_dir, file_name)
        for iris in annotation["iris"]:
            commute_dm.utils.add_annotation_to_file(
                iris,
                annotation["entity_id"],
                input_file_path,
                output_file_path,
            )
    # for dir_ in [COVID_DM_OUTPUT_DIR]:
    #     # for dir_ in [PD_DM_OUTPUT_DIR, COVID_DM_OUTPUT_DIR]:
    #     for file_name, file_path in utils.list_dir(dir_):
    #         for old_label, new_label in label_to_label.items():
    #             utils.replace_label(
    #                 f'name="{old_label}"',
    #                 f'name="{new_label}"',
    #                 file_path,
    #                 file_path,
    #             )
    #             utils.replace_label(
    #                 f">{old_label}<",
    #                 f">{new_label}<",
    #                 file_path,
    #                 file_path,
    #             )
