%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from masstodon.read.txt import spectrum_from_txt
from masstodon.masstodon import masstodon_single, load_masstodon

datapath = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP/wave_height_1.5/FRL-010513-SUBP-WH 1,5-WV 10.txt"
mz, intensity = spectrum_from_txt(datapath)

from masstodon.data.substance_p import substance_p
from masstodon.data.constants import infinity


mz, intensity = substance_p['spectrum']

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3


m = masstodon_single(mz, intensity, fasta, q,
                     min_mz_diff   = infinity,
                     modifications = modifications,
                     orbitrap      = False,
                     threshold     = threshold,
                     isotopic_coverage = isotopic_coverage,
                     min_prob = min_prob, 
                     std_cnt  = std_cnt)

m.dump("dump")


# m.plotly("dump")
# m.plotlygl("dump", shape='rectangles')
# m.good_mols


# m.imperator.solutions[0].plot()
# m.imperator.solutions[0].plot_fittings()
# xx = m.imperator.used_peak_groups()

# m.groups.min_mz[xx] 
# m.groups.max_mz[xx]

prec = list(m.ome.sources())[0]
for mol in m.ome.observables():
    name = m.ome.G[prec][mol]['name']
    if name == "precursor":
        print(name, mol)


prec_source = next(m.ome.sources())
prec_source.formula
prec_source.modifications

pure_prec = next(m.ome.observables())
iso       = pure_prec.isotopologues()
pure_prec.formula

from masstodon.plot.spectrum import plot_spectrum
import matplotlib.pyplot as plt


local_total_intensity = m.spec[1342:1360].total_intensity()
max_intensity = max(m.spec[1342:1360].intensity)
colors = ["red", "blue", "green"]

m.spec.plot(show=False)
i = 0
for mol in m.ome.observables():
    name = m.ome.G[prec][mol]['name']
    if name == "precursor" and mol.q == 1:
        iso = mol.isotopologues()
        plt.scatter(x = iso.mz,
                    y =iso.probability * .5 * local_total_intensity,
                   color=colors[i])
        i += 1
plt.show()


# plot bands of tolerance around the experimental peaks.
plt.vlines(m.groups.min_mz, 0, max_intensity, colors='red', linestyle='dotted')
plt.vlines(m.groups.max_mz, 0, max_intensity, colors='blue', linestyle='dotted')
i = 0
for mol in m.ome.observables():
    name = m.ome.G[prec][mol]['name']
    if name == "precursor" and mol.q == 1:
        iso = mol.isotopologues()
        plt.scatter(x = iso.mz,
                    y =iso.probability * .5 * local_total_intensity,
                   color=colors[i])
        i += 1
m.spec.plot()
# all seems to fall in correct positions.


for o in m.ome.observables():
    print(m.ome.G[prec][o]['name'], o, o.isotopologues().mz)





pure_prec.monoisotopic_mz

iso.mz
iso.probability
m.mz_digits


subspec = m.spec[1346:1360]
subspec.mz


m.threshold
m.spec.threshold

m.groups.min_mz
m.groups.max_mz

# testing the bloody plotly bar


