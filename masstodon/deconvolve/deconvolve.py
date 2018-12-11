from    networkx.linalg.attrmatrix  import attr_matrix
import  numpy                       as     np
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass

from masstodon.models.nnls import nnls
from masstodon.models.quantile_deconvolution import quant_deconv

# TODO: sort out the interpretation of betas problem.
class DeconvolutionProblem(object):
    def fit(self, 
            connected_component,
            total_intensities,
            min_mz,
            max_mz,
            mean_mz,
            include_zero_intensities = False):
        self.include_000 = include_zero_intensities
        self.cc     = connected_component
        X, ordering = attr_matrix(self.cc, edge_attr='prob')
        ordering    = np.array(ordering)
        mol_columns = ordering < 0
        peak_rows   = ordering >= 0
        X           = X[:,mol_columns][peak_rows,:]
        self.idx    = ordering[peak_rows]
        self.total_intensities = total_intensities[self.idx]
        self.mol_idx = ordering[mol_columns]
        self.mz_s   = min_mz[self.idx]
        self.mz_e   = max_mz[self.idx]
        self.mean_mz= mean_mz[self.idx]
        Y           = self.total_intensities
        if self.include_000:
            Y1 = np.concatenate((Y, np.zeros(X.shape[1])))
            x  = 1.0 - np.array(X.sum(axis=0)).flatten()
            X1 = np.concatenate((X, np.diag(x)))
            self.model = self.fit_model(X1, Y1)
        else:
            self.model = self.fit_model(X, Y)
        self.X = X
        self.Y = Y

    def fit_model(self, X, Y):
        raise NotImplementedError

    def XY(self):
        """Return contingency matrix and response vector."""
        return self.X, self.Y

    def iter_estimates(self):
        yield from zip(self.mol_idx, self.model.coef())

    def __len__(self):
        return len(self.cc)

    def spectrum(self):
        pred = self.model.fitted()
        Y = self.Y
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

    def l1_abs(self):
        return self.model.l1_abs()

    def l1_rel(self):
        return self.model.l1_rel()

    def total_intensity(self):
        return sum(self.Y)

    def fitted(self):
        # TODO: wtf?
        # if self.include_000:
        #     fitted = self.model.fitted()
        #     return fitted
        # else:
        #     return self.model.fitted()
        return self.model.fitted()

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

    def plot_fittings(self,
                      plt_style = 'fast',
                      peaks     = True,
                      show      = True):
        plt.style.use(plt_style)
        x     = self.mean_mz
        y_hat = self.model.fitted()
        if peaks:
            plt.vlines(x, [0], y_hat, color='red')
        plt.scatter(x, y_hat, s = 16, color='red')
        if show:
            plt.show()


class NNLSDeconvolution(DeconvolutionProblem):
    def fit_model(self, X, Y):
        return nnls(X, Y)

    def __repr__(self):
        return "Ich bin eine Nonnegetive Regression. Hast du Angst?"


class QuantileDeconvolution(DeconvolutionProblem):
    def fit_model(self, X, Y):
        return quant_deconv(X, Y)

    def __repr__(self):
        return "Ich bin eine Quantil Regression. Hast du Angst?"



def deconvolve(connected_component,
               total_intensities,
               min_mz,
               max_mz,
               mean_mz,
               method = 'nnls',
               include_zero_intensities = False):
    if method == 'nnls':
        dp = NNLSDeconvolution()
    elif method == 'quantile':
        dp = QuantileDeconvolution()
    else:
        raise NotImplementedError("The is not such method of deconvolution as {}.".format(method))

    dp.fit(connected_component,
           total_intensities,
           min_mz,
           max_mz,
           mean_mz,
           include_zero_intensities)
    return dp

