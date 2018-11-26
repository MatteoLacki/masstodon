%load_ext autoreload
%autoreload 2

from masstodon.precursor.precursor import precursor

P = precursor('RPKPQQFFGLM', 3, name='substanceP')
P[3]
P.c_fragments()
P.molecules()

mols = list(P.molecules())
mols[40][0].formula

list(P.c_fragments())
list(P.z_fragments())
list(P.uncharged_molecules())

P.name