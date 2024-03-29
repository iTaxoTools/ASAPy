from itaxotools.common.param import Field, Group

def params():
	return Group(key='root', children=[
		Group(key='general', label='General', children=[
			Field(key='all',
				label='Generate all files',
				doc=("Generate all partition files."),
				type=bool,
				default=True
			),
			Field(
				key = 'mega',
				label = 'MEGA CSV',
				doc = ("Set if the provided distance Matrix is CSV file\nexported from MEGA 6 or X."),
				type = bool,
				default = False
			),
			Field(
				key = 'sequence_length',
				label = 'Sequence length',
				doc = ("Original length of the sequence."),
				type = int,
				default = 600
			),
		]),
		Group(key='advanced', label='Advanced', children=[
			Field(key='number',
				label='Scores Kept',
				doc=("Number of results with the highest scores to be displayed."),
				type=int,
				default=10
			),
			Field(key='seuil_pvalue',
				label='Probability',
				doc=("Limit for results to be reported."),
				type=float,
				default=0.01
			),
			Field(key='seed',
				label='Seed',
				doc=("Use fixed seed value. If you don’t want to use a fixed seed value, set to -1."),
				type=int,
				default=-1
			),
		]),
		Group(key='distance', label='Distance', children=[
			Field(key='method',
				label='Method',
				doc=("If you provide a fasta file, you can select the substitution model to compute the distances."),
				type=int,
				list={
                    0: 'Kimura-2P (K80)',
                    1: 'Jukes-Cantor (JC69)',
                    2: 'Tamura-Nei (TN93)',
                    3: 'Simple Distance',
				},
				default=1
			),
			Field(key='rate',
				label='Kimura TS/TV',
				doc=("Transition/transversion for Kimura 3-P distance."),
				type=float,
				default=2.0
			),
		]),
	])

# 
# params = {
#     "general": {
#         "label": "General",
#         "fields": {
#             "all": {
#                 "label":    "Generate all files",
#                 "doc":      "Generate all partition files.",
#                 "type":     "bool",
#                 "default":  True
#             },
#             "mega": {
#                 "label":    "MEGA CSV",
#                 "doc":      "Set if the provided distance Matrix is CSV file\nexported from MEGA 6 or X.",
#                 "type":     "bool",
#                 "default":  False
#             },
#             "sequence_length": {
#                 "label":    "Sequence length",
#                 "doc":      "Original length of the sequence.",
#                 "type":     "int",
#                 "default":  600
#             },
#         }
#     },
#     "advanced": {
#         "label": "Advanced",
#         "doc": "Advanced",
#         "fields": {
#             "number": {
#                 "label":    "Scores Kept",
#                 "doc":      "Number of results with the highest scores to be displayed.",
#                 "type":     "int",
#                 "default":  10
#             },
#             "seuil_pvalue": {
#                 "label":    "Probability",
#                 "doc":      "Limit for results to be reported.",
#                 "type":     "float",
#                 "default":  0.01
#             },
#             "seed": {
#                 "label":    "Seed",
#                 "doc":      "Use fixed seed value. If you don’t want to use a fixed seed value, set to -1.",
#                 "type":     "int",
#                 "default":  -1
#             },
#             # "replicates": {
#             #     "label":    "Replicates",
#             #     "doc":      "Number of replicates for statistical tests.",
#             #     "type":     "int",
#             #     "default":  1000
#             # },
#             # "pond_pente": {
#             #     "label":    "Pond Pente",
#             #     "doc":      "Internal parameter, for testing only.",
#             #     "type":     "float",
#             #     "default":  0.1
#             # },
#         }
#     },
#     "distance": {
#         "label": "Distance",
#         "fields": {
#             "method": {
#                 "label":    "Method",
#                 "doc":      "If you provide a fasta file, you can select the substitution model to compute the distances.",
#                 "type":     "list",
#                 "default":  1,
#                 "data": {
#                 "items":  [0,1,2,3],
#                 "labels": ["Kimura-2P (K80)","Jukes-Cantor (JC69)","Tamura-Nei (TN93)","Simple Distance"]
#                 }
#             },
#             "rate": {
#                 "label":    "Kimura TS/TV",
#                 "doc":      "Transition/transversion for Kimura 3-P distance.",
#                 "type":     "float",
#                 "default":  2.0
#             }
#         }
#     }
# }
