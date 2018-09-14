"""Ome is more than proteome, or genome, or metalobome. It's simply ome."""
import intervaltree         as iTree
import networkx             as nx
import matplotlib.pyplot    as plt
from collections import defaultdict

from masstodon.precursor.precursor  import precursor, FalsePrecursor
from masstodon.molecule.molecule    import molecule
from masstodon.formula.formula      import Formula, dict2string


class Ome(object):
    def __init__(self, iso_calc):
        self.iso_calc = iso_calc

    def make_graph(self, precursors=[], molecules=[]):
        self.G = nx.Graph()
        for p_kwds in precursors:
            prec = precursor(iso_calc = self.iso_calc, **p_kwds)
            self.G.add_node(prec, type='source')
            for mol, name in prec.molecules():
                self.G.add_node(mol, type='observable')
                self.G.add_edge(mol, prec, name=name)
        for m_kwds in molecules:
            name = m_kwds['name']
            m_kwds = m_kwds.copy()
            del m_kwds['name']
            mol = molecule(iso_calc=self.iso_calc, **m_kwds)
            prec = FalsePrecursor(name=name,
                                  iso_calc=self.iso_calc,
                                  **m_kwds)
            self.G.add_node(prec, type='source')
            self.G.add_node(mol,  type='observable')
            self.G.add_edge(mol, prec, name='precursor')

    def plot(self, plt_style = 'seaborn-poster',
                   node_label_size = 10,
                   node_size = 10,
                   show      = True):
        plt.style.use(plt_style)
        layout = nx.spring_layout(self.G)
        node_labels = defaultdict(str)
        for s in self.sources():
            node_labels[s] = s.name
        node_labels = nx.draw_networkx_labels(self.G, 
                                              pos       = layout,
                                              labels    = node_labels,
                                              font_size = node_label_size)
        nx.draw_networkx_nodes(self.G, pos=layout, node_size = node_size)
        nx.draw_networkx_edges(self.G, pos=layout)
        edge_labels = {(a,b): self.G[a][b]['name'] for a,b in self.G.edges}
        nx.draw_networkx_edge_labels(self.G,
                                     pos=layout,
                                     edge_labels=edge_labels)
        if show:
            plt.show()

    def iter_nodes(self, that_type):
        for n in self.G:
            if self.G.node[n]['type'] is that_type:
                yield n

    def observables(self):
        yield from self.iter_nodes('observable')

    def sources(self):
        yield from self.iter_nodes('source')

    def write(self, path):
        pass

    def filter_by_deviations(self, 
                             subspectra,
                             std_cnt = 3):
        """Filter good subspectra and molecules.

        Good are those, that have a chance to be close to each other.
        This is estimated by the N-sigma rule: 3 sigmas should approximately
        contain 99.9 percent of probability.
        """
        emp_tree = iTree.IntervalTree()
        for subspec in subspectra:
            s, e = subspec.interval
            emp_tree[s:e] = subspec
        good_mols = []
        good_subspectra = set([])
        for mol in self.observables():
            s, e = mol.interval(std_cnt = std_cnt)
            touched_spectra = emp_tree[s:e]
            if touched_spectra:
                status = 'ok'
                good_mols.append(mol)
                good_subspectra |= touched_spectra
            else:
                status = 'std_dev'
            self.G.node[mol]['status'] = status
        return good_mols, good_subspectra
