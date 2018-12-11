%load_ext autoreload
%autoreload 2

import pandas   as pd
from os         import listdir, makedirs
from os.path    import join as pjoin, exists as pexists

from masstodon.read.txt         import spectrum_from_txt
from masstodon.masstodon        import masstodon_single
from masstodon.data.constants   import infinity

common_path = "/Users/matteo/Projects/masstodon/data/Belgian/2013_05_01/SUBP"
folders = ("wave_velocity_300", "wave_height_1.5")

def iter_data(path, folders):
    """Iterate over existing spectra in different folders with common root at the path."""
    for f in folders:
        pf = pjoin(path, f)
        for cff in listdir(pf):
            cffp = pjoin(pf, cff)
            try:
                yield cff, spectrum_from_txt(cffp)
            except Exception as e:
                pass


dump_folder = "/Users/matteo/Projects/masstodon/dumps/belgian/synapt"

threshold         = 0.075
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
timings           = True
include_zero_intensities = False

def get_WH_WV(f):
    f = f.replace(" ","")
    f = f.split(".")[0]
    WH, WV = f.split("-")[3:5]
    WH = WH[2:]
    if WH == '1,5':
        WH = 1.5
    else:
        WH = int(WH)/100
    WV = int(WV[2:])
    return WH, WV

bad = []
molecules = []
# i, (cff, (mz, intensity)) = next(enumerate(iter_data(common_path, folders)))
def filter_by_WV(wv_value):
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        wh, wv = get_WH_WV(cff)
        if wv == wv_value:
            yield i, (cff, (mz, intensity))
            break
# i, (cff, (mz, intensity)) = next(filter_by_WV(400))

def res_with_figs():
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        try:
            m, t = masstodon_single(mz, intensity, fasta, q,
                                 min_mz_diff        = infinity,
                                 modifications      = modifications,
                                 orbitrap           = False,
                                 threshold          = threshold,
                                 isotopic_coverage  = isotopic_coverage,
                                 min_prob           = min_prob, 
                                 std_cnt            = std_cnt,
                                 get_timings        = timings,
                                 include_zero_intensities = include_zero_intensities)
            cff = cff.replace(".txt", "")
            df = pjoin(dump_folder, "{i}_".format(i=i) + cff)
            if not pexists(df):
                makedirs(df)
            m.dump(df)
            m.write(df)
            m.plotly(pjoin(df, "spectrum.html"),shape='rectangles',show=False)
            # WWW = enumerate(m.ome.iter_molecule_estimates())
            # j,M = next(WWW)
            # j,M = next(WWW)
            wh, wv = get_WH_WV(cff)
            for j, M in enumerate(m.ome.iter_molecule_estimates()):
                if j > 0:
                    molecules.append((M[0],M[1],M[2],M[3],M[4],M[5][0],cff,wh,wv,i))
            print(f'Finished with {cff}')
        except AssertionError:
            bad.append((cff, mz, intensity))

res_with_figs()
D = pd.DataFrame(molecules)
D.columns = ('formula','q','g','qg_formula','intensity',
             'name','experiment','WH','WV','index')
D.to_csv(pjoin(dump_folder,"intensities.csv"), index=False)

