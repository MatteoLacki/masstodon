"""A class to create subproblems for deconvolution.

Should be in principle replaced by a simple C++ code.
Now, it really should."""

import  networkx            as      nx
from    collections         import  defaultdict, Counter
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass


from masstodon.deconvolve.deconvolve import deconvolve 


class Imperator(object):
    def __init__(self, molecules,
                       clustering,
                       min_prob = 0.8,
                       isotopic_coverage = .99):
        """Long live the Empire."""
        for prob in (min_prob, isotopic_coverage):
            assert prob > 0.0 and prob <= 1.0
        self.min_prob = min_prob
        self.P        = isotopic_coverage
        self.clust    = clustering
        self.mols     = molecules

    def divide(self):
        self.G = nx.Graph()
        for tree in self.divide_iter():
            self.G.add_edges_from((M_cnt, E, {'prob': P}) for 
                                  (M_cnt, E), P in tree.items())

    def divide_iter(self):
        ls = self.clust.ls
        for M_cnt, M in enumerate(self.mols):
            M_cnt = - M_cnt - 1 # it is quicker not to run the __hash__ for M, but use this count
            P_within_groups = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.P, True):
                E = ls[I_mz]
                if E > 0:
                    I[(M_cnt, E)]   += I_prob
                    P_within_groups += I_prob
            if P_within_groups >= self.min_prob: 
                # planting tree E - M - E in graph G
                yield I

    def impera_iter(self):
        """Iterate over connected components of the deconvolution graph."""
        return nx.connected_component_subgraphs(self.G)

    def impera(self):
        """List all conected components."""
        self.solutions = [deconvolve(cc, 
                                     self.clust.groups.intensity,
                                     self.clust.groups.min_mz,
                                     self.clust.groups.max_mz,
                                     self.clust.groups.mean_mz) for cc in self.impera_iter()]

    def set_estimated_intensities(self):
        for sol in self.solutions:
            for idx, estimate in sol.iter_estimates():
                if estimate > 0.0:
                    idx = - idx- 1
                    self.mols[idx].intensity = estimate

    def plot(self,
             node_size = 2.0,
             show      = True):
        """Plot the deconvolution graph."""
        colors = ['red' if N > 0 else 'blue' for N in self.G]
        # red  : peak groups
        # blue : theoretical molecules  
        nx.draw(self.G,
                node_size  = node_size,
                node_color = colors)
        if show:
            plt.show()

    def plot_ccs(self,
                 plt_style = 'dark_background',
                 show      = True):
        """Plot sizes of the small deconvolution graphs."""
        plt.style.use(plt_style)
        edges_cnt = []
        mols_cnt  = []
        clust_cnt = []
        for sol in self.solutions:
            edges_cnt.append(len(sol.cc.edges))
            x = Counter(0 if N < 0 else 1 for N in sol.cc)
            mols_cnt.append(x[0])
            clust_cnt.append(x[1])
        plt.scatter(mols_cnt, clust_cnt, s = edges_cnt)
        plt.title("Analysis of Connected Components")
        plt.xlabel("Number of Molecules")
        plt.ylabel("Number of Peak Groups")
        plt.axes().set_aspect('equal')
        if show:
            plt.show()

    def plot_solutions(self, fancy     = True,
                             plt_style = 'fast',
                             bar_color = 'grey',
                             bar_alpha = 0.5,
                             show      = True):
        plt.style.use(plt_style)
        for sol in self.solutions:
            plotter = sol.plot_fancy if fancy else sol.plot
            plotter(plt_style, bar_color, bar_alpha, False)
        if show:
            plt.show()

    def __len__(self):
        return len(self.G)


def divide_ed_impera(molecules,
                     clustering,
                     min_prob          = .8,
                     isotopic_coverage = .99):
    imp = Imperator(molecules, clustering, min_prob, isotopic_coverage)
    imp.divide()
    imp.impera()
    imp.set_estimated_intensities()
    return imp


