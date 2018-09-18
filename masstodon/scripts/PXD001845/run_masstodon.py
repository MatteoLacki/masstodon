%load_ext autoreload
%autoreload 2
%load_ext line_profiler


from itertools  import islice
import os
from os.path    import join as pjoin, exists as pexists
from time       import time

from masstodon.scripts.PXD001845.iter_folders import path
from masstodon.scripts.PXD001845.iter_folders import get_charge
from masstodon.scripts.PXD001845.iter_folders import iter_scans
from masstodon.masstodon import masstodon_single


path_iter = iter_scans(path)

fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
min_prob          = .8
isotopic_coverage = .999
std_cnt           = 3
dump_path         = "/Users/matteo/Projects/masstodon/dumps/one_process"


# the simplest case: no multiple threads yet.
# mz, intensity, q, path = next(path_iter)
for mz, intensity, q, path in islice(path_iter, 1):
    try:
        t0 = time()
        M = masstodon_single(mz, intensity, fasta, q, '', 
                             isotopic_coverage = isotopic_coverage,
                             min_prob          = min_prob, 
                             std_cnt           = std_cnt)
        local_dump_path = pjoin(dump_path,
                                pjoin(*path.split('/')[-2:]))
        if not pexists(local_dump_path):
            os.makedirs(local_dump_path)
        M.dump(local_dump_path)
        M.write(local_dump_path)
        t1 = time()
    except Exception as e:
        print(path)
        print(e)

M.ome.G_stats


M.dump("dump")
M.imperator.errors_to_json("dump/errors.json")



import numpy as np

dp = M.imperator.solutions[1]
dp.model.l1()
dp.model.l2()
dp.model.l1_error()
dp.model.l2_error()
dp.model.l1()/dp.model.Y.sum()
np.sum(np.abs(dp.model.res()))
np.sum(np.abs(dp.model.res()))



dp.model.l1()

dp.l2()
dp.Y
dp.Y.sum()
dp.plot()
dp.l1_error()

import numpy as np

used = set([])
used.update(p for s in M.imperator.solutions for p in s.idx)
used = np.array(list(used))
used.sort()

group_intensity = M.spec.bc.groups.intensity


used = np.zeros(shape=group_intensity.shape, dtype=bool)
used[used_idx] = True

group_intensity[used].sum() / group_intensity[~used].sum()



len(M.spec.bc.groups.intensity)


M.imperator.solutions[0].idx

M.spec.plot(show=False)
M.imperator.plot_solutions()
M.spec.l1()
M.spec.l2()

sum(dp.l1() for dp in M.imperator.solutions)/M.spec.l1()
sum(dp.l2() for dp in M.imperator.solutions)/M.spec.l2()

M.ome.G_stats

[n for n in M.ome.G]
[e for e in M.ome.G.edges]

