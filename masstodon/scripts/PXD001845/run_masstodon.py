%load_ext autoreload
%autoreload 2
%load_ext line_profiler


from itertools  import islice
import json
import os
from os.path    import join as pjoin, exists as pexists
from time       import time

from masstodon.scripts.PXD001845.iter_folders import path
from masstodon.scripts.PXD001845.iter_folders import get_charge
from masstodon.scripts.PXD001845.iter_folders import iter_scans, non_modified_scans
from masstodon.masstodon import masstodon_single

# def iter_data(path):
#     for mz, intensity, charge, experiment in iter_scans(path):
#         exp = experiment.split("/")[-2]
#         if "AMB_Bora" in exp and "precZ" in exp and "ETD" in exp:
            # yield mz, intensity, charge, experiment
def iter_data(path):
    for mz, intensity, charge, experiment in iter_scans(path):
        exp = experiment.split("/")[-2]
        if "AMB_pBora":
            yield mz, intensity, charge, experiment


if __name__ == "__main__":
    # path_iter = iter_scans(path)
    data = list(iter_data(path))
    fasta             = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
    min_prob          = .8
    isotopic_coverage = .999
    std_cnt           = 3
    dump_path         = "/Users/matteo/Projects/masstodon/dumps/one_process"
    stop     = None
    max_iter = len(data) if stop is None else min(stop, len(data))
    i = 0
    # the simplest case: no multiple threads yet.
    # mz, intensity, q, path = next(path_iter)
    for mz, intensity, q, path in islice(data, stop):
        i += 1
        print(f"Spec {i} out of {max_iter}.")
        try:
            M, timings = masstodon_single(mz, intensity, fasta, q, '', 
                                          isotopic_coverage = isotopic_coverage,
                                          min_prob          = min_prob, 
                                          std_cnt           = std_cnt,
                                          get_timings       = True)
            local_dump_path = pjoin(dump_path,
                                    pjoin(*path.split('/')[-2:]))
            if not pexists(local_dump_path):
                os.makedirs(local_dump_path)
            M.dump(local_dump_path, indent=4)
            M.write(local_dump_path)
            M.plotlygl(local_dump_path, show=False)
            with open(pjoin(local_dump_path, 'timings.json'), 'w') as h:
                json.dump(timings, h, indent=4)
            print("\tOK\n")
        except Exception as e:
            print(path)
            print(e)
