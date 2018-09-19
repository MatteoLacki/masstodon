"""A class to create subproblems for deconvolution.

Should be in principle replaced by a simple C++ code.
Now, it really should."""
try:
    import plotly
    import plotly.graph_objs as go
    plotly_available = True
except ImportError:
    plotly_available = False
import json
import networkx            as      nx
import numpy               as      np
from   collections         import  defaultdict, Counter
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
        self.isotopic_coverage = isotopic_coverage
        self.clust    = clustering
        self.mols     = molecules

    def divide(self):
        self.G = nx.Graph()
        for tree in self.divide_iter():
            self.G.add_edges_from((M_cnt, E, {'prob': P}) for 
                                  (M_cnt, E), P in tree.items())

    def save_graph(self, path):
        """Save the computed deconvolution graph.

        Do assure that the folders in the 'path' do exist beforehand.

        Parameters
        ==========
        path : str
            Path where to store the deconvolution graph.
        """
        nx.write_gpickle(self.G, path)

    def load_graph(self, path):
        """Load a computed graph.

        Parameters
        ==========
        path : str
            Path where the deconvolution graph is stored.
        """
        self.G = nx.read_gpickle(path)

    def divide_iter(self):
        ls = self.clust.ls
        for M_cnt, M in enumerate(self.mols):
            M_cnt = - M_cnt - 1 # this counter is quicker than the __hash__ for M
            P_within_groups = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.isotopic_coverage, True):
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
        # this set will contain all the indices of used mean_mzs.
        self.used_idx = set([])
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

    def plot_solutions(self, plt_style = 'fast',
                             bar_color = 'grey',
                             bar_alpha = 0.5,
                             show      = True):
        plt.style.use(plt_style)
        h   = self.clust.groups.intensity
        x_l = self.clust.groups.min_mz
        x_r = self.clust.groups.max_mz
        plt.bar(x       = x_l,
                height  = h,
                bottom  = [0],
                width   = x_r - x_l,
                align   = 'edge',
                alpha   = bar_alpha,
                color   = bar_color)
        for sol in self.solutions:
            sol.plot_fittings(plt_style, False)
        if show:
            plt.show()

    def plot_solutions_simple(self, plt_style = 'fast',
                              bar_color = 'grey',
                              bar_alpha = 0.5,
                              show      = True):
        plt.style.use(plt_style)
        h = self.clust.groups.intensity
        x = self.clust.groups.mean_mz
        plt.vlines(x, [0], h, color='grey')
        for sol in self.solutions:
            sol.plot_fittings(plt_style, False, False)
        if show:
            plt.show()

    def plotly_solutions(self, path):
        """Make a plotly plot of the fittings.

        The plot overlays scatterplot of fitted intensities
        over the bars corresponding to total intensities in peak groups.

        Parameters
        ==========
        path : str
            Where to save the plot.
            The file should have 'html' extension to avoid warnings.
        """
        # some simplifications
        if plotly_available:
            # data 4 plot
            x = self.clust.groups.mean_mz
            h = self.clust.groups.intensity.astype(int)
            l = self.clust.groups.min_mz
            r = self.clust.groups.max_mz
            fitted_int = []
            fitted_mzs = []
            fitted_to_int = []
            for s in self.solutions:
                fitted_mzs.extend(list(s.mean_mz))
                fitted_int.extend(list(s.model.fitted()))
                fitted_to_int.extend(list(s.model.Y))
            fitted_int      = np.array(fitted_int)
            fitted_mzs      = np.array(fitted_mzs).astype(int)
            fitted_to_int   = np.array(fitted_to_int).astype(int)
            text_annotation = np.array([f"fit {fit_int:.0f}<br>obs {fit2int:.0f}"
                           for fit_int, fit2int in zip(fitted_int, fitted_to_int)])
            # plot elements
            spectrum_bars = go.Bar(x=x, y=h, width=r-l, name="Peak Group")
            fitted_dots = go.Scatter(x = fitted_mzs,
                                     y = fitted_int,
                                     text      = text_annotation, # maps to labels
                                     hoverinfo = "text",          # show only labels
                                     mode      = "markers",       # default='lines'
                                     marker    = {"color" : "orange"},
                                     name      = "Fitted")
            data = [spectrum_bars, fitted_dots]
            layout = go.Layout(
                title = "Observed versus Fitted Spectrum",
                font  =dict(
                    color="white"
                ),
                yaxis = dict(
                    title = "Intensity",
                    color = "white"
                ),
                xaxis = dict(
                    title = "mass/charge [Th]",
                    color = "white"
                ),
                plot_bgcolor = "black",
                paper_bgcolor= "black"
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename=path)
        else:
            raise ImportError("You must install plotly seperately! We want you, the 'enlighted coder', to have the possibility to run masstodon with pypy out-of-the-box.")


    def __len__(self):
        return len(self.G)

    def solutions_l1_error_abs(self):
        """Summarize the l1 error for all the solutions"""
        return sum(s.model.l1_abs() for s in self.solutions)

    def solutions_total_intensity(self):
        return sum(s.model.total_intensity() for s in self.solutions)

    def solutions_l1_error_rel(self):
        return self.solutions_l1_error_abs() /\
              (self.total_fitted() + self.total_intensity_fitted_to())

    def total_intensity(self):
        return sum(self.clust.groups.intensity)

    def used_peak_groups(self):
        """Get a mask selecting peak groups used in the fitting."""
        used_idx = set()
        used_idx = set([])
        used_idx.update(p for s in self.solutions for p in s.idx)
        used_idx = np.array(list(used_idx))
        used = np.zeros(shape = self.clust.groups.intensity.shape,
                        dtype = bool)
        used[used_idx] = True
        return used

    def total_intensity_fitted_to(self):
        return sum(self.clust.groups.intensity[self.used_peak_groups()])

    def total_fitted(self):
        return sum(s.model.total_fitted() for s in self.solutions)

    def l1_abs(self):
        intensity = self.clust.groups.intensity
        intensity_beyond_theory = sum(intensity[~self.used_peak_groups()])
        return intensity_beyond_theory + self.solutions_l1_error_abs()

    def l1_rel(self):
        return self.l1_abs()/(self.total_intensity() + self.total_fitted())

    def errors_to_json(self, path):
        errors = {}
        errors['l1_abs'] = self.l1_abs()
        errors['l1_rel'] = self.l1_rel()
        errors['solutions_l1_error_abs'] = self.solutions_l1_error_abs()
        errors['solutions_l1_error_rel'] = self.solutions_l1_error_rel()
        with open(path, 'w') as f:
            json.dump(errors, f)


def imperator(molecules,
              clustering,
              min_prob          = .8,
              isotopic_coverage = .99):
    imp = Imperator(molecules, clustering, min_prob, isotopic_coverage)
    imp.divide()
    imp.impera()
    imp.set_estimated_intensities()
    return imp


def load_imperator(molecules,
                   clustering,
                   deconvolution_graph_path,
                   min_prob          = .8,
                   isotopic_coverage = .99):
    imp = Imperator(molecules, clustering, min_prob, isotopic_coverage)
    imp.load_graph(deconvolution_graph_path)
    imp.impera()
    imp.set_estimated_intensities()
    return imp