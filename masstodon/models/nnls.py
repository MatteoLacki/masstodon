try:
    import matplotlib.pyplot as plt
except RuntimeError:
    print('Warning: Python not installed as framework. You will not be able to use matplotlib.')
    pass
import numpy as np
import scipy.optimize

from masstodon.models.model import Model


class NNLS(Model):
    def fit(self, X, Y):
        self.X = X
        self.Y = Y
        self._coef, self._res = scipy.optimize.nnls(X, Y)

    def l2(self):
        return self._res

    def coef(self):
        return self._coef

    def __call__(self, x):
        return np.dot(x, self._coef)

    def predict(self, x):
        return self(x)

    def res(self):
        return self.Y - self.fitted()

    def __repr__(self):
        return "Ich bin eine Nichtnegative Regression. Hast du Angst?"

    def fitted(self):
        return np.array(np.dot(self.X, self._coef)).flatten()

    def plot(self, 
             plt_style = 'default',
             lines_col = 'grey',
             show      = True):
        plt.style.use(plt_style)
        x = range(len(self.Y))
        fitted = self.fitted()
        plt.vlines(x = x, ymin=[0], ymax=self.Y, colors=lines_col)
        plt.scatter(x, fitted, color = 'blue')
        plt.scatter(x, self.Y, color = 'red')
        plt.title("Analysis of fitting:\nRed: real Y, Blue: fitted values.")
        plt.xlabel("Index")
        plt.ylabel("Y values")
        if show:
            plt.show()

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

def nnls(X, Y):
    model = NNLS()
    model.fit(X, Y)
    return model
