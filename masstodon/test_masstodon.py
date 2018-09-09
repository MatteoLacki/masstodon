# from collections import defaultdict, Counter
# import unittest

# from masstodon.data.get_dataset import get_dataset
# from masstodon.masstodon import masstodon


# class TestMassTodonOnSubstanceP(unittest.TestCase):
#     def setUp(self):
#         """Load data on Substance P."""

#         substanceP = get_dataset('substanceP')
#         modifications = defaultdict(dict)
#         for (number, group), mods in substanceP.precursor.modifications.items():
#             modifications[number][group] = dict(mods)
        
#         self.args = {'spectrum':        substanceP.spectrum,
#                      'min_intensity':   50.0,
#                      'fasta':           substanceP.precursor.fasta,
#                      'name':            'substanceP',
#                      'modifications':   modifications,
#                      'charge':          3,
#                      'mz_tol':          .05}
        
#     def run_masstodon(self):
#         """Run one session of the MassTodon on the Substance P data."""
#         masstodon = MassTodon(**self.args)

#     def test_some_times(self):
#         """Test the MassTodon on substance P, M times."""
#         M = 100
#         res = Counter() 
#         for i in range(M):
#             try:
#                 masstodon = self.run_masstodon()
#                 res['success'] += 1
#             except ValueError as e:
#                 res['failure'] += 1
        
#         self.assertEqual(res['success'], M)
#         self.assertEqual(res['failure'], 0)
                                                  
