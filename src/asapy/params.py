params = {
    "general": {
        "label": "General",
        "fields": {
            "mega": {
                "label":    "MEGA CVS",
                "doc":      "If True, the distance Matrix is supposed to be MEGA CVS\n(other formats are guessed).",
                "type":     "bool",
                "default":  False
            },
            "all": {
                "label":    "Generate all files",
                "doc":      "Generate all partition and tree files.",
                "type":     "bool",
                "default":  True
            },
            "sequence_length": {
                "label":    "Sequence length",
                "doc":      "Original length of the sequence.",
                "type":     "int",
                "default":  600
            },
        }
    },
    "advanced": {
        "label": "Advanced",
        "doc": "Advanced",
        "fields": {
            "replicates": {
                "label":    "Replicates",
                "doc":      "Number of replicates for statistical tests.",
                "type":     "int",
                "default":  1000
            },
            "seuil_pvalue": {
                "label":    "Pvalue Threshold",
                "doc":      "Limit for results to be reported.",
                "type":     "float",
                "default":  0.05
            },
            "pond_pente": {
                "label":    "Pond Pente",
                "doc":      "Limit for results to be reported.",
                "type":     "float",
                "default":  0.1
            },
            "seed": {
                "label":    "Seed",
                "doc":      "Use fixed seed value. If you don’t want to use a fixed seed value, set to -1.",
                "type":     "int",
                "default":  -1
            },
        }
    },
    "distance": {
        "label": "Distance",
        "fields": {
            "method": {
                "label":    "Method",
                "doc":      "If you provide a fasta file, you can select the substitution model to compute the distances.",
                "type":     "list",
                "default":  1,
                "data": {
                "items":  [0,1,2,3],
                "labels": ["Kimura-2P (K80)","Jukes-Cantor (JC69)","Tamura-Nei (TN93)","Simple Distance"]
                }
            },
            "rate": {
                "label":    "Kimura TS/TV",
                "doc":      "Transition/transversion for Kimura 3-P distance.",
                "type":     "float",
                "default":  2.0
            }
        }
    }
}
