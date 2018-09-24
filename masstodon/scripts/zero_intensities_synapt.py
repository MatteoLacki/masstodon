%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import pandas   as pd
from   os       import listdir, makedirs
from   os.path  import join as pjoin, exists as pexists

from masstodon.read.txt       import spectrum_from_txt
from masstodon.masstodon      import masstodon_single
from masstodon.data.constants import infinity

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

dump_folder = "/Users/matteo/Projects/masstodon/dumps/belgian/synapt"

threshold         = 0.025
fasta             = 'RPKPQQFFGLM'
modifications     = {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}}
q                 = 3
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
timings           = True
include_zero_intensities = True

bad = []


# i, (cff, (mz, intensity)) = next(enumerate(iter_data(common_path, folders)))
def iter_outcomes(include_zero_intensities):
    for i, (cff, (mz, intensity)) in enumerate(iter_data(common_path, folders)):
        try:
            WH, WV = get_WH_WV(cff)
            m, t = masstodon_single(mz, intensity, fasta, q,
                                 min_mz_diff        = 1.1,
                                 modifications      = modifications,
                                 orbitrap           = False,
                                 threshold          = threshold,
                                 isotopic_coverage  = isotopic_coverage,
                                 min_prob           = min_prob, 
                                 std_cnt            = std_cnt,
                                 get_timings        = timings,
                                 include_zero_intensities = include_zero_intensities)
            df = pjoin(dump_folder, str(i) + "_" + cff)
            if not pexists(df):
                makedirs(df)
            m.dump(df)
            m.write(df)
            m.plotlygl(df, shape='rectangles', show=False)
            row = {"i":i, "exp": cff, "WH": WH, "WV": WV}
            row.update(m.imperator.errors())
            row.update({"t_" + str(n): T for n,T in t})
            row.update(m.ome.G_stats)
            for s in ('ETDorHTR', 'ETnoD_PTR_fragments', 'ETnoD_precursor', 'PTR_precursor'):
                row["cz_simple." + str(s)] = int(m.cz_simple.intensities[s])
                row["cz.{s}" + str(s)]        = int(m.cz.intensities[s])
            for i in range(11):
                row["cz_simple.prob." + str(i)] = m.cz_simple.probabilities["fragmentation_bond"].get(i, 0.0)
                row["cz.prob." + str(i)] = m.cz.probabilities["fragmentation_bond"].get(i, 0.0)
            yield row
            print('Finished with '+ str(cff))
        except AssertionError:
            bad.append((cff, mz, intensity))

results_for_plot = pd.DataFrame(iter_outcomes(True))
results_for_plot.to_csv(
    "/Users/matteo/Projects/masstodon/dumps/belgian/synapt_results.csv",
    index = False)
results_for_plot_with000 = results_for_plot.copy()


results_for_plot = pd.DataFrame(iter_outcomes(False))
results_for_plot.to_csv(
    "/Users/matteo/Projects/masstodon/dumps/belgian/synapt_results.csv",
    index = False)
results_for_plot_ohne000 = results_for_plot.copy()

pd.DataFrame({'ohne':results_for_plot_ohne000['nodes_nonzero_intensity'],
              'with':results_for_plot_with000['nodes_nonzero_intensity']})
