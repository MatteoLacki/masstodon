ptms = {
"phosphorylation":              {"P":1, "H":1,  "O":3},
"acetylation":                  {"C":2, "H":2,  "O":1},
"glycation_by_hexose":          {"C":6, "H":10, "O":5},
"amidation":                    {"H":1, "N":1,  "O":-1},
"hydroxylation":                {"O":1},
"methylation":                  {"C":1, "H":2},
"ubiquitylation":               {"C":378, "H":627, "N":105, "O":117, "S":1},
"pyroglutamic_acid":            {"H":-2, "O":-1},
"sulfation":                    {"S":1, "O":3},
"gamma_carboxyglutamic_acid":   {"C":1, "O":2},
"dehydralation":                {"H":-2, "O": -1},
"oxidation":                    {"O":1},
"carboxymethylation":           {"C":2, "H":2, "O":2},
"palmitic_acylation":           {"C":16,"H":30,"O":1},
"oleic_acylation":              {"C":18,"H":32,"O":1},
"arachidonic_acylation":        {"C":20,"H":30,"O":1},
"docosahexanoic_acylation":     {"C":22,"H":30,"O":1}
}

typical_ptm_sites = {
    "phosphorylation": ("S", "T", "Y")
}