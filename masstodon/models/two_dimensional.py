import numpy as np
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass

from masstodon.models.model import Model


class Model2D(Model):
    def x_percentiles(self, No):
        """Get percentiles of 'x'.
        
        Parameters
        ----------
        No : int
            Order of percentiles, e.g. No = 4 corresponds to quartiles, 5 to quintiles, 10 to deciles.
        """
        return np.percentile(self.x, np.linspace(0, 100, No+1))

    def plot(self,
             knots_no  = 1000,
             plot_data = True,
             plt_style = 'dark_background',
             scatter_color = 'blue',
             show      = True):
        """Plot data points and the fitted line.

        Parameters
        ----------
            knots_no : int
                Number of inner knots for the line.
            plot_data : logical
                Include the fit-to-data in the plot.
            plt_style : str
                The type of the matplotlib style used in the plot.
            colors_no : int
                The number of colors to cycle through.
            show : logical
                Immediately show the plot? Alternatively, just add it to the present canvas.

        """
        plt.style.use(plt_style)
        if plot_data:
            plt.scatter(self.x, self.y, c=scatter_color, s=.5)
        x_ = np.linspace(self.x_min, self.x_max, knots_no)
        y_ = self(x_)
        plt.plot(x_, y_, c='red')
        if show:
            plt.show()

    def res(self):
        """Get residuals/errors of the model."""
        return self.y - self(self.x)

    def error_stats(self):
        err         = self.res()
        err_deciles = np.percentile(err, np.linspace(0, 100, 11))

        median_absolute_deviation = np.median(np.abs(err))
        return err_deciles, median_absolute_deviation

