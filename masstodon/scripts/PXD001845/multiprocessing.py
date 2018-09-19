import multiprocessing
import gc
import time

no_cores = 3

tasks   = multiprocessing.Queue()
returns = multiprocessing.Queue()


def subworker(q):
    q.get()
    gc.collect()
    gc.disable()
    t   = time.time()
    ret = do_work()
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
        p.join(5)
            # Block the calling thread until the process whose join() 
            # method is called terminates or until the OPTIONAL TIMEOUT occurs.

        if p.is_alive():
            p.terminate()
            returns.put((task, 'timeout'))
        else:
            returns.put((task, q.get()))


procs = [ multiprocessing.Process(target=worker) for _ignored in range(no_cores)]

# so you can put task after preparing the workers.

for x in range(20):
    tasks.put(x)


for proc in procs:
    proc.start()

for x in range(no_cores):
    tasks.put('end')


for proc in procs:
    proc.join()


while not returns.empty():
    print(returns.get())

