from masstodon.scripts.PXD001845.csv2fasta import csv2fasta

class Suspicion(Exception):
    pass


def read_fasta(folder):
    matching_csv = [f for f in csv2fasta if folder in f]
    if len(matching_csv) == 1:
        fasta = csv2fasta[matching_csv[0]]
        return fasta
    else:
        raise Suspicion("Suspicious to have more than one folder here.")
