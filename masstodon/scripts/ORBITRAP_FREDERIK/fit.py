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

class WrongSystem(Exception):
    pass
if platform == "darwin": # latte on skimmed soya milk anyone?
    common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/UBI_ORBI/12plus precursor"
    dump_folder = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/res/12plus"
    out_json = "/Users/matteo/Projects/masstodon/data/Belgian/2015_07/res/12plus/ubi.json"
    processes_no = 4
elif platform == "linux": # death metal and long dirty hair, fuck yeah!
    common_path = "/mnt/disk/masstodon/data/belgian/2015_07/UBI_ORBI/12plus"
    dump_folder = "/mnt/disk/masstodon/data/belgian/2015_07/UBI_ORBI/12plusRes"
    out_json = "/mnt/disk/masstodon/data/belgian/2015_07/UBI_ORBI/12plusRes.json"
    processes_no = 24
else:
    raise WrongSystem("Path not specified correctly.")

files = [f for f in listdir(common_path) if ".mzXML" in f]

ubiquitin='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'

def iter_data(max_files_no=None, max_scans_no=None):
    for f_idx, f in enumerate(islice(files, max_files_no)):
        for x in islice(read_mzml(pjoin(common_path,f)), max_scans_no):
            if not x['centroided'] and x['msLevel']==2:
                scan = int(x['num'])
                mz = x['m/z array']
                intensity = x['intensity array']
                yield f_idx, f, scan, mz, intensity

# f_idx, f, scan, mz, intensity = next(iter_data(1,1))
def single_run(f_idx, f, scan, mz, intensity):
    f = f.replace(".mzXML","")
    param = float(f.split("_")[-1].replace(",","."))
    row = {"f_idx":f_idx, "exp": f, "scan":scan, "param": param}
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
        df = pjoin(dump_folder, str(f_idx)+"_"+f, str(scan))
        if not pexists(df):
            makedirs(df)
        m.dump(df)
        m.write(df)
        m.plotly(pjoin(df, "spectrum.html"), show=False)
        row.update(m.imperator.errors())
        row.update({"t_"+str(n): T for n,T in t})
        row.update(m.ome.G_stats)
        row['estimates'] = list(m.ome.iter_molecule_estimates())
        row['success'] = True
        print('Finished with {}'.format(f))
    except Exception as e:
        row['success'] = False
        print("Error in file {} in scan {}".format(f, scan))
    return row

# data_it = iter_data(2,2)
data_it = iter_data()

T0 = time()
with mp.Pool(processes_no) as p:
    stats = p.starmap(single_run, data_it)
T1 = time()
with open(out_json, "w") as f:
    json.dump(stats, f, indent=4)
print("Total fit time for all ETD spectra: {t}".format(t=T1-T0))
