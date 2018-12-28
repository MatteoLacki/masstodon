%load_ext autoreload
%autoreload 2

from itertools import islice
from os import listdir, makedirs
from os.path import join as pjoin, exists as pexists
from sys import platform
import json
import multiprocessing as mp
from time import time

from masstodon.data.constants import infinity
from masstodon.masstodon import masstodon_single
from masstodon.read.mzml import read_mzml
from masstodon.plot.spectrum import plot_spectrum

common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/UBI_ORBI/12plus precursor"
dump_folder = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/res/12plus"
out_json = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/res/12plus/ubi.json"
processes_no = 4

f = "FRL_220715_ubi_714_ETD_0,03.mzXML"
scan = 133
ubiquitin='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

def iter_data(f, max_scans_no=None):
    for x in islice(read_mzml(pjoin(common_path,f)), max_scans_no):
        if not x['centroided'] and x['msLevel']==2:
            scan = int(x['num'])
            mz = x['m/z array']
            intensity = x['intensity array']
            yield f, scan, mz, intensity

f, scan, mz, intensity = next(filter(lambda x: x[1] == scan,
                                     iter_data(f)))
# f, scan, mz, intensity = next(iter_data(f))
# plot_spectrum(mz, intensity)

def single_run(f, scan, mz, intensity):
    f = f.replace(".mzXML","")
    param = float(f.split("_")[-1].replace(",","."))
    row = {"scan":scan, "param": param}
    df = pjoin(dump_folder, str(0)+"_"+f, str(scan))
    try:
        print('Running {}'.format(f))
        m, t = masstodon_single(
            mz, intensity, ubiquitin, 12,
            orbitrap           = True,
            threshold          = "0.05 Th",
            isotopic_coverage  = .999,
            min_prob           = .8, 
            std_cnt            = 3,
            get_timings        = True)
        if not pexists(df):
            makedirs(df)
        m.dump(df)
        m.write(df)
        m.plotly(pjoin(df, "spectrum.html"), show=True)
        row.update(m.imperator.errors())
        row.update({"t_"+str(n): T for n,T in t})
        row.update(m.ome.G_stats)
        row['estimates'] = list(m.ome.iter_molecule_estimates())
        row['success'] = True
        print('Finished with {}'.format(f))

        m.precursors

    except Exception as e:
        row['success'] = False
    return row, m

r, m = single_run(f, scan, mz, intensity)

data_it = iter_data(f)

T0 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, data_it)
T1 = time()


len(stats)
stats[0]['success']

from collections import Counter
Counter(s['success'] for s in stats)

errors = list(filter(lambda x: not x['success'], stats))


from masstodon.precursor.precursor  import precursor
prec = precursor(fasta=ubiquitin,
                 q=12,
                 iso_calc=m.iso_calc)
it = prec.molecules()
prec_mol = list(filter(lambda x: x[1] == 'precursor' and x[0].q == 12, prec.molecules()))[0]
env = prec_mol[0].isotopologues()

m.spec

import matplotlib.pyplot as plt
plot_spectrum(mz, intensity, show=False)
plt.scatter(env.mz, env.probability*sum(intensity)/4)
plt.show()


