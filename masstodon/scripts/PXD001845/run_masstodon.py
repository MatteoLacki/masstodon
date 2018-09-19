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


stop = None
# the simplest case: no multiple threads yet.
# mz, intensity, q, path = next(path_iter)
for mz, intensity, q, path in islice(path_iter, stop):
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
        M.plotly(local_dump_path, show=False)
        t1 = time()
    except Exception as e:
        print(path)
        print(e)

