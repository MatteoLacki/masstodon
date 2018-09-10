# """Testing the setting up of the deconvolution problems."""
# from   collections import Counter
# import networkx as nx
# import unittest

# from ..data.get_dataset  import get_dataset
# from .deconvolve         import deconvolve
# from ..spectrum.spectrum import spectrum


# class TestDeconvolution(unittest.TestCase):
#     def test_glue_sister_isotopologues(self):
#         """Test if glueing isotopologues works."""
#         G = nx.Graph()
#         G.add_node('M0')
#         G.add_node('I0', mz=1.1, probability=.1); G.add_edge('M0','I0')
#         G.add_node('I1', mz=1.2, probability=.05); G.add_edge('M0','I1')
#         G.add_node('G0'); G.add_edge('I0','G0'); G.add_edge('I1','G0')
#         G.add_node('I2', mz=2.1, probability=.2); G.add_edge('M0','I2')
#         G.add_node('I3', mz=2.2, probability=.25); G.add_edge('M0','I3')
#         G.add_node('I4', mz=2.3, probability=.15); G.add_edge('M0','I4')
#         G.add_node('G1'); G.add_edge('I2','G1'); G.add_edge('I3','G1'); G.add_edge('I4','G1')
#         # GRAPH LOOKS LIKE THIS
#         # M-\----\--\--\
#         # |  \    \  \  \
#         # I0 I1   I2 I3 I4
#         #  \ /      \ | /
#         #  G0        G1
#         _glue_sister_isotopologues(G)
#         for g in G:
#             if g[0] == 'G':
#                 for I in G[g]:
#                     self.assertEqual(len(G[I]), 2)
#                     for M in G[I]:
#                         if M[0] == 'M':
#                             self.assertEqual(M, 'M0')
#         E_stats = dict(Counter(N[0] for N in G))
#         R_stats = {'G':2, 'I':2, 'M':1}
#         self.assertEqual(E_stats, R_stats)

#     def test_deconvolve(self):
#         """Test the numbers of nodes in different connected componets."""
#         # R_ = real
#         # nodes:    M  I   G    M  I   G    M  I   G
#         R_stats = {(1, 7, 7), (2, 14, 8), (3, 21, 9)}   # _merge_sister_Is=True
#         # R_stats = {(1, 11, 7), (2, 27, 8), (3, 54, 9)}# _merge_sister_Is=False

#         subP = get_dataset('substanceP')
#         # There are 6 precursors
#         precursors = list(m for m in subP.precursor.molecules()
#                           if m.name == 'precursor')

#         spectrum = sum(m.isotopologues() for m in precursors)
#         spectrum = Spectrum(mz=spectrum.mz,
#                             intensity=100000 * spectrum.probability)
#         spectrum.round_mz(precision=2)
#         DGs = deconvolve(precursors,
#                          spectrum,
#                          'Matteo',
#                          mz_tol=.05,
#                          isospec_args={'mz_digits': 2},
#                          _merge_sister_Is=True)
#         # E_ = expected
#         E_stats = [Counter(N[0] for N in DG) for DG in DGs ]
#         E_stats = set([(s['M'], s['I'], s['G']) for s in E_stats])
#         self.assertEqual(R_stats, E_stats)

# if __name__ == "__main__":
#     unittest.main()
