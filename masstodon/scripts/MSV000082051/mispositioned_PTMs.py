from masstodon.precursor.precursor import precursor
from masstodon.data.amino_acids import amino_acids

from masstodon.isotopes.calculator import isotope_calculator
from masstodon.formula.formula import NegativeAtomCount

amino_acids[0]
iso = isotope_calculator()

try:
    prec = precursor(fasta="A",
                     q=3,
                     iso_calc=iso,
                     name='test',
                     modifications={1:{"C_alpha":{"O":-2}}})
except NegativeAtomCount:
    print("Got you!!!")

for p_kwds in experiments[X]['precursors']:
    try:
        prec = precursor(fasta="A",
                         q=3,
                         iso_calc=iso,
                         name='test',
                         modifications={1:{"C_alpha":{"O":-2}}})
    except NegativeAtomCount:
        print("Got you!!!")

it = filter_PTM_assignments([
    dict(fasta="A",
         q=3,
         iso_calc=iso,
         name='test',
         modifications={1:{"C_alpha":{"O":-2}}}),
    dict(fasta="Q",
         q=3,
         iso_calc=iso,
         name='test',
         modifications={1:{"C_alpha":{"O":-1}}})
])
list(it)

def filter_ptm_assignments(precursors_dictionaries):
    for prec_kwds in precursors_dictionaries:
        try:
            prec = precursor(**prec_kwds)
            yield prec_kwds
        except NegativeAtomCount:
            pass
