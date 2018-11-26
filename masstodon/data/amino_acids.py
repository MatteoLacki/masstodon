from masstodon.formula.formula import Formula


#TODO: simplify this for the sake of fuck!
amino_acids = [['A',
                [['C_alpha', [['H', 4], ['C', 2]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['C',
                [['C_alpha', [['H', 4], ['C', 2], ['S', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['E',
                [['C_alpha', [['H', 6], ['C', 4], ['O', 2]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['D',
                [['C_alpha', [['H', 4], ['C', 3], ['O', 2]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['G',
                [['C_alpha', [['H', 2], ['C', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['F',
                [['C_alpha', [['H', 8], ['C', 8]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['I',
                [['C_alpha', [['H', 10], ['C', 5]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['H',
                [['C_alpha', [['H', 6], ['C', 5], ['N', 2]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['K',
                [['C_alpha', [['H', 11], ['C', 5], ['N', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['M',
                [['C_alpha', [['H', 8], ['C', 4], ['S', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['L',
                [['C_alpha', [['H', 10], ['C', 5]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['N',
                [['C_alpha', [['H', 5], ['C', 3], ['O', 1], ['N', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['Q',
                [['C_alpha', [['H', 7], ['C', 4], ['O', 1], ['N', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['P',
                [['C_alpha', []],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 7], ['C', 4], ['N', 1]]]]],
               ['S',
                [['C_alpha', [['H', 4], ['C', 2], ['O', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['R',
                [['C_alpha', [['H', 11], ['C', 5], ['N', 3]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['T',
                [['C_alpha', [['H', 6], ['C', 3], ['O', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['W',
                [['C_alpha', [['H', 9], ['C', 10], ['N', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['V',
                [['C_alpha', [['H', 8], ['C', 4]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]],
               ['Y',
                [['C_alpha', [['H', 8], ['C', 8], ['O', 1]]],
                 ['C_carbo', [['C', 1], ['O', 1]]],
                    ['N', [['H', 1], ['N', 1]]]]]]


def get_amino_acids():
    """Retrieve the information on amino acidic bricks.

    Returns
    =======
    A dictionary with keys (element, backbone_atom_group).
    The values are linear counters storing atom counts.
    Possible backbone_atom_group include: N, C_carbo, C_alpha.
    """
    AAS = {(aa_name, brick): Formula({a: c for a, c in atom_cnt})
           for aa_name, bricks in amino_acids
           for brick, atom_cnt in bricks}
    return AAS

