# from collections import Counter, namedtuple
# import unittest

# from masstodon.data.get_dataset             import get_dataset
# from masstodon.estimates_matcher.cz_simple  import SimpleCzMatch
# from masstodon.estimates_matcher.cz         import CzMatch
# from masstodon.reporter.misc                import Brick

# # TODO add test for the optimality of the fragment matching

# def node_to_tuple(node):
#     return node.type + str(node.no), node.q

# class VoidClass(object):
#     pass


# no alphas now: rewrite test using graphs
#class TestPeakPicker(unittest.TestCase):
#    def setUp(self):
#        """Set up a method."""
#
#        mol = get_dataset('substanceP')
#        mols = list(mol.precursor.molecules())
#        mols_dict = {(m.name, m.q, m.g): m for m in mols}
#        results = [
#            {'alphas': [
#                {'estimate': 1000,
#                 'molecule': mols_dict[('c4', 1, 0)]},
#                {'estimate': 5000,
#                 'molecule': mols_dict[('c5', 1, 0)]},
#                {'estimate': 8000,
#                 'molecule': mols_dict[('c6', 1, 0)]}],
#             'status': 'optimal'},
#            {'alphas': [
#                {'estimate': 1200,
#                 'molecule': mols_dict[('z7', 1, 0)]},
#                {'estimate': 4800,
#                 'molecule': mols_dict[('z6', 1, 0)]},
#                {'estimate': 8300,
#                 'molecule': mols_dict[('z5', 1, 0)]}],
#             'status': 'optimal'},
#            {'alphas': [
#                {'estimate': 40000,
#                 'molecule': mols_dict[('precursor', 3, 0)]},
#                {'estimate': 30000,
#                 'molecule': mols_dict[('precursor', 2, 0)]},
#                {'estimate': 20000,
#                 'molecule': mols_dict[('precursor', 1, 0)]}],
#             'status': 'optimal'}]
#
#        molecules = [m['molecule'] for res in results for m in res['alphas']]
#        # for res in results:
#        #     for x in res['alphas']:
#        #         mol = VoidClass()
#        #         mol.intensity = x['estimate']
#        #         x['molecule'].intensity =
#        #         brick = Brick(molecule=x['molecule'])
#        #         bricks.append(brick)
#        #
#        # # Preparing a phoney MassTodon
#        # self.precursor = VoidClass()
#        self.precursor_charge = 3
#        self.molecules = molecules
#
#    def tearDown(self):
#        """Tear down a method."""
#        pass
#
#    def test_SimpleCzMatch(self):
#        print("Testing SimpleCzMatch graph creation.")
#
#        # precursor_charge = 3
#        matches = SimpleCzMatch(molecules=self.molecules,
#                                precursor_charge=self.precursor_charge)
#
#        expected_nodes = set([('c4', 1), ('z7', 1), ('c5', 1),
#                              ('z6', 1), ('c6', 1), ('z5', 1)])
#        obtained_nodes = set(node_to_tuple(n)
#                             for n in matches.graph.nodes)
#        self.assertEqual(expected_nodes, obtained_nodes)
#
#        expected_edges = set(frozenset(e) for e in
#                     [[('c4', 1), ('z7', 1)],
#                      [('c5', 1), ('z6', 1)],
#                      [('c6', 1), ('z5', 1)],
#                      [('c4', 1), ('c4', 1)],
#                      [('c5', 1), ('c5', 1)],
#                      [('c6', 1), ('c6', 1)],
#                      [('z5', 1), ('z5', 1)],
#                      [('z6', 1), ('z6', 1)],
#                      [('z7', 1), ('z7', 1)]])
#
#        obtained_edges = set(frozenset([node_to_tuple(n),
#                                        node_to_tuple(m)])
#                             for n, m in matches.graph.edges)
#        self.assertEqual(expected_edges, obtained_edges)
#
    # def test_intermediate_matchmaker(self):
    #     print("Testing basic matchmaker: graph creation.")
    #     matches = czMatchMakerIntermediate(masstodon=self.masstodon)

    #     expected_nodes = set([('c4', 1, 0), ('z7', 1, 0), ('c5', 1, 0),
    #                           ('z6', 1, 0), ('c6', 1, 0), ('z5', 1, 0)])
    #     obtained_nodes = set(matches.graph.nodes)
    #     self.assertEqual(expected_nodes, obtained_nodes)

    #     expected_edges = set(frozenset(e) for e in
    #                         [[('c4', 1, 0), ('z7', 1, 0)],
    #                          [('c5', 1, 0), ('z6', 1, 0)],
    #                          [('c6', 1, 0), ('z5', 1, 0)]])
    #     obtained_edges = set(frozenset(e) for e in matches.graph.edges)
    #     self.assertEqual(expected_edges, obtained_edges)

    # def test_advanced_matchmaker(self):
    #     print("Testing advanced matchmaker: graph creation.")
    #     pass

if __name__ == "__main__":
    unittest.main()
