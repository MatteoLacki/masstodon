"""Ome is more than proteome, or genome, or metalobome. It's simply ome."""
import intervaltree         as iTree
import json
import networkx             as nx
import matplotlib.pyplot    as plt
from math        import ceil
from collections import defaultdict

from masstodon.precursor.precursor  import precursor, FalsePrecursor
from masstodon.molecule.molecule    import molecule
from masstodon.formula.formula      import Formula, dict2string
from masstodon.write.csv_tsv        import write_rows

class Ome(object):
    def __init__(self, iso_calc):
        self.iso_calc = iso_calc

    def make_graph(self, precursors=[], molecules=[]):
        self.G = nx.Graph()
        for p_kwds in precursors:
            prec = precursor(iso_calc = self.iso_calc, **p_kwds)
            self.G.add_node(prec, s=True)
            for mol, name in prec.molecules():
                self.G.add_node(mol, s=False)
                self.G.add_edge(mol, prec, name=name)
        for m_kwds in molecules:
            name   = m_kwds['name']
            m_kwds = m_kwds.copy()
            del m_kwds['name']
            mol  = molecule(iso_calc=self.iso_calc, **m_kwds)
            prec = FalsePrecursor(name=name,
                                  iso_calc=self.iso_calc,
                                  **m_kwds)
            self.G.add_node(prec, s=True)
            self.G.add_node(mol,  s=False)
            self.G.add_edge(mol, prec, name='precursor')
        self.G_stats = {'nodes': len(self.G.nodes),
                        'edges': len(self.G.edges)}

    def dump(self, path):
        """Save the computed sources-observables graph.

        Do assure that the folders in the 'path' do exist beforehand.

        Parameters
        ==========
        path : str
            Path where to store the sources-observables graph.
        """
        nx.write_gpickle(self.G, path)

    def load(self, path):
        """Load a computed graph.

        Parameters
        ==========
        path : str
            Path where the deconvolution graph is stored.
        """
        self.G = nx.read_gpickle(path)

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

    def iter_nodes(self, source=True):
        for n in self.G:
            if self.G.node[n]['s'] is source:
                yield n

    def observables(self):
        yield from self.iter_nodes(False)

    def sources(self):
        yield from self.iter_nodes(True)

    def is_one_precursor(self):
        """Check if there is only one precursor."""
        sources = self.sources()
        prec = next(sources)
        try:
            x = next(sources)
            return False
        except StopIteration:
            return True

    def iter_molecule_estimates(self, header=True):
        """Iterate over molecules with positive estimates."""
        if self.is_one_precursor:
            if header:
                yield ('formula', 'q', 'g',
                       'formula with q and g',
                       'intensity', 'name')
            for m in self.observables():
                yield (str(m.formula), m.q, m.g, 
                       m.formula.str_with_charges(m.q, m.g), 
                       ceil(m.intensity), m.name)
        else:
            if header:
                yield ('formula', 'q', 'g', 
                       'formula with q and g',
                       'intensity')
            for m in self.observables():
                yield (str(m.formula), m.q, m.g,
                       m.formula.str_with_charges(m.q, m.g),
                       ceil(m.intensity))

    def write(self, path):
        """Write precise estimates."""
        write_rows(self.iter_molecule_estimates(), path)

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
        bad_mols = []
        good_subspectra = set([])
        for mol in self.observables():
            s, e = mol.interval(std_cnt = std_cnt)
            touched_spectra = emp_tree[s:e]
            if touched_spectra:
                status = 'ok'
                good_mols.append(mol)
                good_subspectra |= touched_spectra
            else:
                bad_mols.append(mol)
        self.G.remove_nodes_from(bad_mols)
        self.G_stats['nodes_std_dev_filter'] = len(self.G.nodes)
        self.G_stats['edges_std_dev_filter'] = len(self.G.edges)
        return good_mols, good_subspectra

    def filter_by_estimated_intensity(self):
        good = []
        bad  = []
        for o in self.observables():
            if o.intensity > 0.0:
                good.append(o)
            else:
                bad.append(o)
        self.G.remove_nodes_from(bad)
        self.G_stats['nodes_nonzero_intensity'] = len(self.G.nodes)
        self.G_stats['edges_nonzero_intensity'] = len(self.G.edges)

    def dump_stats(self, path):
        with open(path, 'w') as f:
            json.dump(self.G_stats, f)




def ome(iso_calc, precursors=[], molecules=[]):
    """Generate an ome. """
    _ome = Ome(iso_calc)
    _ome.make_graph(precursors, molecules)
    return _ome


def load_ome(iso_calc, path):
    """Load a dumped ome. """
    _ome = Ome(iso_calc)
    _ome.load(path=path)
    return _ome