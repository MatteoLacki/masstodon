import intervaltree as iTree


def filter_subspectra_molecules(subspectra,
                                molecules,
                                std_cnt = 3):
    """Filter good subspectra and molecules.

    Good are those, that have a chance to be close to each other.
    This is estimated by the N-sigma rule: 3 sigmas should approximately
    contain 99.9 percent of probability.
    """
    mols = molecules
    emp_tree = iTree.IntervalTree()
    for subspec in subspectra:
        s, e = subspec.interval
        emp_tree[s:e] = subspec
    good_mols = []
    good_subspectra = set([])
    for mol in mols:
        s, e = mol.interval(std_cnt = std_cnt)
        touched_spectra = emp_tree[s:e]
        if touched_spectra:
            good_mols.append(mol)
            good_subspectra |= touched_spectra
    return good_mols, good_subspectra