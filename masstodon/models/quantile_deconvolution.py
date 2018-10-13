import numpy as np
from scipy.optimize import linprog

from masstodon.models.model import Model


class QuantileDeconvolution(Model):
    def fit(self, X, Y, q=.8, l1=.01):
        """Fit the quantile deconvolution.

        Solves the problem:

        Parameters
        ----------
        X : np.array
            The contingency matrix.
        Y : np.array
            The response variables. Shape must be like (K,), not (K,1).
        q : float
            The quantile to compute.
        l1 : float
            The l1 type penalty.
        """
        assert q >= 0.0
        self.X = X
        self.Y = Y
        self.q = q
        self.l1 = l1
        K = self.observations_cnt()
        D = self.controls_cnt()
        # parametrization: x = [s, t, b]
        # s : the negative parts
        # t : the positive parts
        # b : the coefficients of interest
        self.c = np.concatenate(((1-self.q) * np.ones(shape=(K,1)),
                                 self.q     * np.ones(shape=(K,1)),
                                 self.l1    * np.ones(shape=(D,1))))
        self.c.shape = (self.c.shape[0],)
        self.A_ub = np.vstack((
            np.hstack((-np.identity(K),
                        np.zeros(shape=(K,K)),
                        self.X)),
            np.hstack(( np.zeros(shape=(K,K)),
                       -np.identity(K),
                       -self.X)) ))
        self.b_ub = np.concatenate((self.Y, -self.Y))
        self.model = linprog(c=self.c,
                             A_ub=self.A_ub,
                             b_ub=self.b_ub)

    def coef(self):
        """Get the estimated coefficients."""
        return self.model['x'][-self.controls_cnt():]

    def observations_cnt(self):
        """Get the number of observations in the deconvolution problem."""
        return self.Y.shape[0]

    def controls_cnt(self):
        """Get the number of the controls.

        This equals the number of columns in the contingency matrix X."""
        return self.X.shape[1]

    def predict(self, X):
        """Predict responses with the new set of controls."""
        return X @ self.coef()

    def fitted(self):
        """The values fitted to the responses."""
        # @ is matrix multiplication! True story!
        return np.array(self.X @ self.coef()).flatten()

    def res(self):
        """Calculate the residuals."""
        return self.Y - self.fitted()

    def l1_abs(self):
        """Calculate the l1 absolute error."""
        return sum(np.abs(self.res()))

    def __repr__(self):
        return 'This is Quantile regression.'

    def plot_res(self, 
                 plt_style = 'default',
                 lines_col = 'grey',
                 show      = True):
        plt.style.use(plt_style)
        x   = range(len(self.Y))
        res = self.res()
        plt.vlines(x = x, ymin=[0], ymax=res, colors=lines_col)
        plt.scatter(x, res, color = 'orange')
        plt.title("Analysis of residuals")
        plt.xlabel("Index")
        plt.ylabel("Residual")
        if show:
            plt.show()

def quant_deconv(X, Y, q=.8, l1=.1):
    model = QuantileDeconvolution()
    model.fit(X, Y, q, l1)
    return model
