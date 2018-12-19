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
            if 'q' not in prec_kwds:
                prec_kwds['q'] = 1
            prec = precursor(iso_calc=iso,
                             **prec_kwds)
            yield prec_kwds
        except NegativeAtomCount:
            pass
