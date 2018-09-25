from masstodon.data.ptms import ptms


def modify_fasta(fasta):
    modifications = {pos + 1: {"C_alpha": ptms["phosphorylation"].copy()}
                        for pos, aa in enumerate(fasta) 
                        if aa in ("B", "J")}
    replacements = {'B': 'T', 'J': 'S'}
    for aa in replacements:
        fasta = fasta.replace(aa, replacements[aa])
    return fasta, modifications
