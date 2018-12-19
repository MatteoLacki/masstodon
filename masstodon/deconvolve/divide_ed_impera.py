"""A class to create subproblems for deconvolution.

Should be in principle replaced by a simple C++ code.
Now, it really should."""
try:
    import plotly
    import plotly.graph_objs as go
    from plotly.io import write_json
    from masstodon.plot.plotly import get_black_layout
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


def rectangles(l, r, h):
    X = []
    Y = []
    for a,b,c in zip(l, r, h):
        X.append(a); Y.append(0)
        X.append(a); Y.append(c)
        X.append(b); Y.append(c)
        X.append(b); Y.append(0)
    return np.array(X), np.array(Y)


def triangles(l, r, h):
    X = []
    Y = []
    for a,b,c in zip(l, r, h):
        X.append(a); Y.append(0)
        X.append((a+b)/2.0); Y.append(c)
        X.append(b); Y.append(0)
    X, Y = np.array(X), np.array(Y)
    i = np.argsort(X)
    return X[i], Y[i]


class Imperator(object):
    def __init__(self, molecules,
                       groups,
                       lightweight_spectrum,
                       min_prob = 0.8,
                       isotopic_coverage = .99):
        """Long live the Empire."""
        for prob in (min_prob, isotopic_coverage):
            assert prob > 0.0 and prob <= 1.0
        self.min_prob = min_prob
        self.isotopic_coverage = isotopic_coverage
        self.groups   = groups
        self.ls       = lightweight_spectrum
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
        for M_cnt, M in enumerate(self.mols):
            M_cnt = - M_cnt - 1 # this counter is quicker than the __hash__ for M
            P_within_groups = 0.0
            # edge is an collection of merged isotopologues
            I = defaultdict(float) # values correspond to total probability on that edge. 
            for I_mz, I_prob in M.isotopologues(self.isotopic_coverage, True):
                E = self.ls[I_mz]
                if E > 0:
                    I[(M_cnt, E)]   += I_prob
                    P_within_groups += I_prob
            if P_within_groups >= self.min_prob: 
                # planting tree E - M - E in graph G
                yield I

    def impera_iter(self):
        """Iterate over connected components of the deconvolution graph."""
        return nx.connected_component_subgraphs(self.G)

    def impera(self,
               deconvolution_method = 'nnls',
               include_zero_intensities=False):
        """List all conected components."""
        # this set will contain all the indices of used mean_mzs.
        self.used_idx = set([])
        self.solutions = [deconvolve(cc, 
                                     self.groups.intensity,
                                     self.groups.min_mz,
                                     self.groups.max_mz,
                                     self.groups.mean_mz,
                                     deconvolution_method,
                                     include_zero_intensities) for cc in self.impera_iter()]

    def set_estimated_intensities(self):
        for sol in self.solutions:
            for idx, estimate in sol.iter_estimates():
                if estimate > 0.0:
                    idx = -idx - 1
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
        h   = self.groups.intensity
        x_l = self.groups.min_mz
        x_r = self.groups.max_mz
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
        h = self.groups.intensity
        x = self.groups.mean_mz
        plt.vlines(x, [0], h, color='grey')
        for sol in self.solutions:
            sol.plot_fittings(plt_style, False, False)
        if show:
            plt.show()

    def __get_xhlr(self):
        x = self.groups.mean_mz
        h = self.groups.intensity.astype(int)
        l = self.groups.min_mz
        r = self.groups.max_mz
        return x, h, l, r

    def __get_fitted(self):
        fitted_int = []
        fitted_mzs = []
        fitted_to_int = []
        for s in self.solutions:
            fitted_mzs.extend(list(s.mean_mz))
            fitted_int.extend(list(s.model.fitted()))
            fitted_to_int.extend(list(s.Y))
        fitted_int      = np.array(fitted_int).astype(int)
        fitted_mzs      = np.array(fitted_mzs)
        fitted_to_int   = np.array(fitted_to_int).astype(int)
        text_annotation = np.array(["fit {fit_int:.0f}<br>obs {fit2int:.0f}".format(fit_int=fit_int,
                                                                                    fit2int=fit2int)
                                    for fit_int, fit2int in zip(fitted_int, fitted_to_int)])
        return fitted_int, fitted_mzs, fitted_to_int, text_annotation

    def plotly(self, path='', shape='triangles', webgl=True, show=True):
        """Make a plotly plot of the fittings.

        The plot overlays scatterplot of fitted intensities
        over the bars corresponding to total intensities in peak groups.

        Parameters
        ==========
        path : str
            Where to save the plot.
            The file should have 'html' extension to avoid warnings.
        webgl : boolean
            Should we use WebGL?
        show : boolean
            Show the plot in browser immediately?
        """
        if plotly_available:
            # data 4 plot
            x, h, l, r = self.__get_xhlr()
            shape = rectangles if shape == 'rectangles' else triangles
            X, Y = shape(l, r, h)
            fitted_int, fitted_mzs, fitted_to_int, text_annotation = self.__get_fitted()
            if webgl:
                scatter = go.Scattergl
            else:
                scatter = go.Scatter
            spectrum_bars = scatter(x = X, y = Y, name = "Peak Group")
            fitted_dots = scatter(x         = fitted_mzs,
                                  y         = fitted_int,
                                  text      = text_annotation, # maps to labels
                                  hoverinfo = "text",          # show only labels
                                  mode      = "markers",       # default='lines'
                                  marker    = {"color" : "orange"},
                                  name      = "Fitted")
            data = [spectrum_bars, fitted_dots]
            layout = get_black_layout()
            self.fig = go.Figure(data=data, layout=layout)
            if path:
                plotly.offline.plot(self.fig,
                                    filename  = path,
                                    auto_open = show)
        else:
            raise ImportError("You must install plotly seperately! We want you, the 'enlighted coder', to have the possibility to run masstodon with pypy out-of-the-box.")

    def dump_fig_to_json(self, path):
        """Dump the plotly figure to json.

        Parameters
        ==========
        path : str
            Where to save json with the plot data and layout.
        """
        write_json(self.fig, file=path)

    def __len__(self):
        return len(self.G)

    def solutions_l1_error_abs(self):
        """Summarize the l1 error for all the solutions"""
        return sum(s.l1_abs() for s in self.solutions)

    def solutions_l1_errors_abs(self):
        """Summarize the l1 error for all the solutions"""
        return np.array([s.l1_abs() for s in self.solutions])

    def solutions_l1_errors_rel(self):
        """Summarize the l1 error for all the solutions"""
        return np.array([s.l1_rel() for s in self.solutions])

    def solutions_total_intensity(self):
        return sum(s.total_intensity() for s in self.solutions)

    def solutions_l1_error_rel(self):
        """Compare the experimental end predicted intensity measures in the areas where traces are found."""
        ITE = self.total_fitted() +\
              self.total_intensity_fitted_to()
        if ITE > 0:
            return self.solutions_l1_error_abs()/ITE
        else:
            # there must have been no experimental intensity at all.
            return 1.0

    def total_intensity(self):
        return sum(self.groups.intensity)

    def used_peak_groups(self):
        """Get a mask selecting peak groups used in the fitting."""
        used_idx = set()
        used_idx = set([])
        used_idx.update(p for s in self.solutions for p in s.idx)
        used_idx = np.array(list(used_idx))
        used = np.zeros(shape = self.groups.intensity.shape,
                        dtype = bool)
        if len(used_idx) > 0:
            used[used_idx] = True
        return used

    def total_intensity_fitted_to(self):
        return sum(self.groups.intensity[self.used_peak_groups()])

    def total_fitted(self):
        return sum(s.model.total_fitted() for s in self.solutions)

    def l1_abs(self):
        intensity = self.groups.intensity
        intensity_beyond_theory = sum(intensity[~self.used_peak_groups()])
        return intensity_beyond_theory + self.solutions_l1_error_abs()

    def l1_rel(self):
        return self.l1_abs()/(self.total_intensity() + self.total_fitted())

    def errors(self):
        errors = {}
        errors['l1_abs'] = self.l1_abs()
        errors['l1_rel'] = self.l1_rel()
        errors['solutions_l1_error_abs'] = self.solutions_l1_error_abs()
        errors['solutions_l1_error_rel'] = self.solutions_l1_error_rel()
        return errors

    def errors_to_json(self, path, indent=None):
        errors = self.errors()
        with open(path, 'w') as f:
            json.dump(errors, f, indent=indent)


def imperator(molecules,
              groups,
              lightweight_spectrum,
              min_prob          = .8,
              isotopic_coverage = .99,
              deconvolution_method = 'nnls',
              include_zero_intensities = False):
    imp = Imperator(molecules,
                    groups,
                    lightweight_spectrum,
                    min_prob,
                    isotopic_coverage)
    imp.divide()
    imp.impera(deconvolution_method,
               include_zero_intensities)
    imp.set_estimated_intensities()
    return imp


def load_imperator(molecules,
                   groups,
                   lightweight_spectrum,
                   deconvolution_graph_path,
                   min_prob                 = .8,
                   isotopic_coverage        = .99,
                   deconvolution_method     = 'nnls',
                   include_zero_intensities = False):
    imp = Imperator(molecules,
                    groups,
                    lightweight_spectrum,
                    min_prob,
                    isotopic_coverage)
    imp.load_graph(deconvolution_graph_path)
    imp.impera(deconvolution_method, include_zero_intensities)
    imp.set_estimated_intensities()
    return imp