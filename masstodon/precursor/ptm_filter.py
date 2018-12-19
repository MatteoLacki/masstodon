from masstodon.formula.formula import NegativeAtomCount
from masstodon.isotopes.calculator import isotope_calculator
from masstodon.precursor.precursor import precursor



def filter_ptm_assignments(precursors_dictionaries):
    """Filter precursors with poorly chosen PTMs.

    Poorly chosen PTMs result in negative atom counts.
    Parameters
    ==========
    precursors_dictionaries : list
        A list of dictionaries with parameters for the precursor function.
    """
    iso = isotope_calculator()
    for prec_kwds in precursors_dictionaries:
        try:
            prec = precursor(iso_calc=iso,
                             q=1,
                             **prec_kwds)
            yield prec_kwds
        except NegativeAtomCount:
            pass
