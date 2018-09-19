ptms = {
"phosphorylation":      {"P":1, "H":1,  "O":3},
"acetylation":          {"C":2, "H":2,  "O":1},
"glycosylation_hexose": {"C":6, "H":10, "O":5},
"amidation":            {"H":1, "N":1,  "O":-1},
"hydroxylation":        {"O":1},
"methylation":          {"C":1, "H":2},
"ubiquitylation":       {"C":378, "H":627, "N":105, "O":117, "S":1},
"pyroglutamic_acid":    {"H":-2, "O":-1},
"sulfation":            {"S":1, "O":3},
"gamma_carboxyglutamic_acid": {"C":1, "O":2}}

typical_ptm_sites = {
    "phosphorylation": ("S", "T", "Y")
}