from masstodon.stats.descriptive import mad_denoising

class Denoiser(object):
    pass

class MADdenoiser(Denoiser):
    def __init__(self, std_cnt = 3):
        self.std_cnt = std_cnt

    def fit(self, x):
        self.is_signal, self.top, self.bottom = mad_denoising(x, self.std_cnt)

    def __call__(self, x):
        return (self.bottom <= x) & (x <= self.top)

