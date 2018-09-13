import gc
import multiprocessing
import time

from masstodon.masstodon import masstodon_base, masstodon_base_load

orbi_data = list(iter_scans(path, paths_only=True))
mz_path, intensity_path, charge, scan_dir = orbi_data[0]

no_cores = 3
timeout  = 2

tasks   = multiprocessing.Queue()
returns = multiprocessing.Queue()
for input_tuple in iter_scans(path, paths_only=True):
    tasks.put(input_tuple)
for x in range(no_cores):
    tasks.put('end')

# Interlude: add the bloody csv output of the results.
todon = masstodon_base_load('dump')
todon.spec.plot()

# 


foo = masstodon_base


def subworker(q):
    foo_kwds = q.get()
    gc.collect()
    gc.disable()
    t0  = time.time()
    res = foo(**foo_kwds)
    t1  = time.time() - t
    q.put(ret)


def worker():
    while True:
        task = tasks.get()
        if task == "end":
            return
        q = multiprocessing.Queue()
        q.put(task)
        p = multiprocessing.Process(target = subworker, args=[q])
        p.start()
        p.join(timeout)
        # Block the calling thread until the process whose join() 
        # method is called terminates or until the OPTIONAL TIMEOUT occurs.

        if p.is_alive():
            p.terminate()
            returns.put((task, 'timeout'))
        else:
            returns.put((task, q.get()))

procs = [ multiprocessing.Process(target=worker) for _ in range(no_cores)]


