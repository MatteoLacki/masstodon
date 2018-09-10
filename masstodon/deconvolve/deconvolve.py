from    networkx.linalg.attrmatrix  import attr_matrix
import  numpy                       as     np
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass

from    masstodon.models.nnls import nnls


class DeconvolutionProblem(object):
    def fit(self, 
            connected_component,
            total_intensities,
            min_mz,
            max_mz,
            mean_mz,
            include_zero_intensities = False):
        self.cc     = connected_component
        mol_columns = np.array([N < 0  for N in self.cc])
        peak_rows   = np.array([N >= 0 for N in self.cc])
        X, ordering = attr_matrix(self.cc, edge_attr='prob')
        X           = X[:,mol_columns][peak_rows,:]
        ordering    = np.array(ordering)
        self.idx    = ordering[ordering >= 0]
        self.total_intensities = total_intensities[self.idx]
        self.mz_s   = min_mz[self.idx]
        self.mz_e   = max_mz[self.idx]
        self.mean_mz= mean_mz[self.idx]
        Y           = self.total_intensities
        if include_zero_intensities:
            Y = np.concatenate((Y, np.zeros(X.shape[1])))
            x = 1.0 - np.array(X.sum(axis=0)).flatten()
            X = np.concatenate((X, np.diag(x)))
        self.X = X
        self.Y = Y
        self.model = nnls(X, Y)

    def XY(self):
        """Return contingency matrix and response vector."""
        return self.X, self.Y

    def iter_estimates(self):
        coefs = np.nditer(self.model.coef())
        for N in self.cc:
            if N < 0:
                yield N, float(next(coefs))

    def __len__(self):
        return len(self.cc)

    def l1(self):
        return sum(np.abs(self.model.res()))

    def l2(self):
        """This is copying what I once complained to Ludi about :)"""
        return self.model.l2()

    def spectrum(self):
        pred = self.model.fitted()
        Y    = self.model.Y
        return self.mz_s, self.mz_e, self.mean_mz, pred, Y[Y > 0]

    def plot(self, plt_style = 'fast',
                   bar_color = 'grey',
                   bar_alpha = 0.5,
                   show      = True):
        plt.style.use(plt_style)
        plt.bar(x       = self.mz_s,
                height  = self.model.fitted(),
                width   = self.mz_e - self.mz_s,
                align   = 'edge',
                alpha   = bar_alpha,
                color   = bar_color)
        Y = self.model.Y
        Y = Y[Y>0]
        plt.vlines(self.mean_mz, [0], Y, linestyles='dotted')
        plt.scatter(self.mean_mz, Y, c = 'red', s = 8)
        if show:
            plt.show()

    def plot_fancy(self,
                   plt_style = 'fast',
                   bar_color = 'grey',
                   bar_alpha = 0.5,
                   show      = True):
        plt.style.use(plt_style)
        x_l   = self.mz_s
        x_r   = self.mz_e
        x     = self.mean_mz
        y_hat = self.model.fitted()
        y = self.model.Y
        y = y[y>0]
        bar_top    = np.maximum(y, y_hat)
        bar_bottom = np.minimum(y, y_hat)
        color      = np.full(y.shape, 'blue')
        color[y < y_hat] = 'red'
        plt.bar(x       = x_l,
                height  = bar_top - bar_bottom,
                bottom  = bar_bottom,
                width   = x_r - x_l,
                align   = 'edge',
                alpha   = bar_alpha,
                color   = color)
        plt.vlines(x, [0], bar_bottom, color='grey')
        plt.scatter(x, y, marker = 'o', s = 16, color='black')
        plt.scatter(x, y_hat, marker = 'x', s = 16, color='black')
        if show:
            plt.show()



def deconvolve(connected_component,
               total_intensities,
               min_mz,
               max_mz,
               mean_mz,
               include_zero_intensities = False):
    dp = DeconvolutionProblem()
    dp.fit(connected_component,
           total_intensities,
           min_mz,
           max_mz,
           mean_mz,
           include_zero_intensities)
    return dp

