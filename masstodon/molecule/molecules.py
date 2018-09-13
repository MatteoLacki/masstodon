"""The possibility to make MassTodon find more things."""

%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.masstodon import MasstodonBase

fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
min_prob = .8
isotopic_coverage = .999
std_cnt = 3

todon = MasstodonBase( mz_digits=4,
                       fasta=fasta,
                       charge=charge,
                       name             = "", 
                       modifications    = {},
                       fragments        = "cz",
                       blocked_fragments= ['c0'],
                       block_prolines   = True,
                       distance_charges = 5.,
                       std_cnt          = 3,
                       isotopic_coverage= .999,
                       min_prob         = .7,
                       deconvolution_graph_path = '' )


todon.set_isotopic_calculator()
todon.set_molecules()
todon.prec
todon.mols

todon.iso_calc

from masstodon.precursor.precursor  import precursor, FalsePrecursor
from masstodon.molecule.molecule    import molecule
from masstodon.formula.formula      import Formula, dict2string

f_dict = {'C':100, 'H':202}
f = Formula(f_dict)
str(f)
f_str = dict2string(f)
f = Formula(f_str)
f.str_with_charges(q=1, g=10)

precursors = [
    dict(fasta = "AAAACCCKKK",
         name  = "a-zine",
         q     = 4),
    dict(fasta = "VTRSQLM",
         q     = 5,
         name  = "b-zine"),
    dict(fasta = "QQQQNNNKKKHIFFDDD",
         q     = 3,
         name  = "c-zine"),
]
molecules = [
    dict(name   = 'saliva',
         formula= 'C100H202',
         q      =  1),
    dict(name   = 'shit',
         formula= 'C102H205',
         q      =  2),
    dict(name   = 'lipstick',
         formula= 'C100H202',
         q      =  2),
    dict(name   = "startrek's leather underwear",
         formula= 'C100H202',
         q      =  2)
]

# OK, so these should be lists of arguments for the respective classes
import intervaltree as iTree
import networkx as nx
import matplotlib.pyplot as plt

class Molecules(object):
    def __init__(self, iso_calc):
        self.iso_calc = iso_calc

    def make_graph(self, precursors, molecules):
        self.G = nx.Graph()
        for p_kwds in precursors:
            prec = precursor(iso_calc = self.iso_calc, **p_kwds)
            self.G.add_node(prec, type='substance')
            for mol in prec.molecules():
                self.G.add_node(mol, type='observable')
                self.G.add_edge(mol, prec)
        for m_kwds in molecules:
            name = m_kwds['name']
            m_kwds = m_kwds.copy()
            del m_kwds['name']
            mol  = molecule(name = 'precursor', iso_calc = self.iso_calc, **m_kwds)
            prec = FalsePrecursor(name = name,  iso_calc = self.iso_calc, **m_kwds)
            self.G.add_node(prec, type='substance')
            self.G.add_node(mol,  type='observable')
            self.G.add_edge(mol, prec)

    def plot(self, plt_style = 'seaborn-poster',
                   node_size = 10,
                   show      = True):
        plt.style.use(plt_style)
        layout = nx.spring_layout(self.G)
        nx.draw_networkx_nodes(self.G, pos=layout, node_size=node_size)
        nx.draw_networkx_edges(self.G, pos=layout)
        if show:
            plt.show()

    def iter_nodes(self, that_type):
        for n in self.G:
            if self.G.node[n]['type'] is that_type:
                yield n

    def observables(self):
        yield from self.iter_nodes('observable')

    def substances(self):
        yield from self.iter_nodes('substance')

    def write(self, path):
        pass

    # def filter_by_deviations(self, 
    #                          subspectra,
    #                          std_cnt = 3):
    #     """Filter good subspectra and molecules.

    #     Good are those, that have a chance to be close to each other.
    #     This is estimated by the N-sigma rule: 3 sigmas should approximately
    #     contain 99.9 percent of probability.
    #     """
    #     mols = molecules
    #     emp_tree = iTree.IntervalTree()
    #     for subspec in subspectra:
    #         s, e = subspec.interval
    #         emp_tree[s:e] = subspec
    #     good_mols = []
    #     good_subspectra = set([])
    #     for mol in mols:
    #         s, e = mol.interval(std_cnt = std_cnt)
    #         touched_spectra = emp_tree[s:e]
    #         if touched_spectra:
    #             good_mols.append(mol)
    #             good_subspectra |= touched_spectra
    #     return good_mols, good_subspectra

mols = Molecules(todon.iso_calc)
mols.make_graph(precursors, molecules)
list(mols.observables())
list(mols.substances())
