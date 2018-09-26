from    math    import  inf
import  numpy   as      np

from masstodon.plot.spectrum import plot_spectrum


def max_x_diff_iterator(x, min_x_diff=1.5):
    """Cluster points based on the distance between them.

    Parameters
    ==========
    x: np.array
        Values to cluster.
    min_x_diff : float
        The minimal m/z difference that separates clusters.

    Yields
    ======
    int : number of cluster
    """
    c  = -1   # guardian
    x_ = -inf # guardian
    for _x in x:
        if _x - x_ >= min_x_diff:
            c += 1
        yield c
        x_ = _x


# around a second. That's something to rewrite later to C++.
def bitonic_iterator(x, w, min_x_diff=.15):
    """Cluster based on bitonicity of the intensity.

    If three consecutive intensities form a dump,
    I_{i-1} > I_i < I_{i+1}, then 'i' is considered
    a new cluster.
    Also, if y_i - y_{i-1} > min_y_diff,
    then 'i' is considered a new cluster.

    Parametes
    ---------
    x: np.array
        Values to cluster, sorted!
    w: np.array
        Weights of y.
    min_x_diff : float
        The minimal difference in x entries that separates clusters.
    """
    # previous two intensities: set up guardians
    i__ = 1.0
    _i_ = 0.0
    c   = -1    # cluster no: the first cluster will get tag '0'
    x_  = -inf # previous m/z 
    for _x, __i in zip(x, w):
        big_x_diff = _x - x_ > min_x_diff
        if (__i >= _i_ and _i_ < i__) or big_x_diff:
            c += 1
        yield c
        if big_x_diff:
            i__ = 0.0
        else:
            i__ = _i_
        _i_ = __i
        x_  = _x


def clusters2array(x, w,
                   clustering=bitonic_iterator,
                   **clustering_kwds):
    """Write elements of cluster iterator to a numpy array.

    Parameters
    ----------
    x: np.array
        Values to cluster, sorted!
    w: np.array
        Weights of y.
    clustering: iterator
        Iterator of clusters.
    *clustering args:
        Unnamed list of arguments for the iterator.
    **clustering_kwds:
        Dictionary of pairs argument-value for the iterator.

    Returns
    -------
        np.array : an array of cluster assignments.
    """
    return np.fromiter(clustering(x, w, **clustering_kwds),
                       dtype=int, count=len(x))


def list_of_clusters(x, w,
                     clustering = bitonic_iterator,
                    *clustering_args,
                   **clustering_kwds):
    """Write elements of cluster iterator to a list.

    Parameters
    ----------
    x: np.array
        Values to cluster, sorted!
    w: np.array
        Weights of y.
    clustering: iterator
        Iterator of clusters.
    *clustering args:
        Unnamed list of arguments for the iterator.
    **clustering_kwds:
        Dictionary of pairs argument-value for the iterator.

    Returns
    -------
        list : list of cluster assignments.
    """
    return list(clustering(x, w, **clustering_kwds))


def iter_cluster_ends(assignments):
    """Iterate over left and right ends of subsequent clusters defined by assignments.

    Provides iteration over 'sparse' representation of the assignments to clusters,
    i.e. tuples composed of the index of the first element in cluster,
    index of the last element of cluster.
    Clusters are assumed to appear orderly one after another.

    Parameters
    ----------
        assignments : iter of int
            Iterator with clusters' assignents.
    Yields
    ------
        tuple: indices marking the beginning and the end of a cluster.
    """
    if isinstance(assignments, np.ndarray):
        c_ = assignments[0]
        assignments = np.nditer(assignments)
        _ = next(assignments)
    else: # an iterable
        c_ = next(assignments)
    i_ = 0
    _i = 1
    for _c in assignments: # next cluster
        if _c == c_ + 1:
            yield i_, _i # start and end of cluster c_
            c_ += 1
            i_ = _i
        _i += 1
    yield i_, _i # start and end of the final cluster


def iter_clusters(x, w, assignments):
    """Iterate over clusters given by assignments.

    Parameters
    ----------
    x: np.array
        Values to cluster, sorted!
    w: np.array
        Weights of y.
    assignments : iter of int
        Iterator with clusters' assignents.
    Yields
    ------
        tuple: indices marking the beginning and the end of a cluster.
    """
    for s, e in iter_cluster_ends(assignments):
        yield x[s:e], w[s:e]


def fix_local_clustering(x, s, e, c,
                         abs_perc_dev = .2):
    """Fix bitonic clustering.

    Assigns the end of a cluster to the next cluster based on 
    the median of distances between peaks.

    Parameters
    ----------
    x :np.array
        The values to be clustered.
    s : int
        The index of the start of the cluster of peaks.
    e :int
        The index of the end of a cluster of peaks.
    c :int
        The number of current cluster of peaks.
    abs_perc_dev : float
        How big m/z deviation is tolerable.
    Yields
    ------
    np.array: assignment into clusters.
    """
    clusters = np.full(shape = (e-s,), fill_value = c, dtype=int)
    x_local = x[s:e]
    # differences of consecutive m/z values
    x_diffs = np.diff(x_local)
    if len(x_diffs) >= 4:
        # the median should be a stable values to compare to
        # as there are vastly more similar diffs than other.
        me_x_diff = np.median(x_diffs)
        cc = c
        signal_border = abs_perc_dev * me_x_diff
        for i, x_diff in enumerate(x_diffs):
            if abs(x_diff - me_x_diff) > signal_border:
                cc += 1
            clusters[i+1] = cc
    return clusters


def bitonic_clustering(x, w,
                       min_x_diff   = .15,
                       abs_perc_dev = .2):
    """Cluster points based on differences in consecutive values of sorted x and the bitonicity of the y.

    Parameters
    ----------
    x: np.array
        Values to cluster, sorted!
    w: np.array
        Weights of y.
    min_x_diff : float
        The minimal difference in x that separates clusters.

    Returns
    -------
    np.array of ints: assignments into consecutive clusters.
    """
    # an array to store the clusters
    clusters = np.full(shape = (len(x),), fill_value = 0, dtype = int)
    # the clustering iterator
    bc = bitonic_iterator(x, w, min_x_diff = min_x_diff)
    # cluster count
    c = 0
    for s, e in iter_cluster_ends(bc):
        lclust = fix_local_clustering(x, s, e, c, abs_perc_dev)
        clusters[s:e] = lclust
        if lclust[0] != lclust[-1]:
            c = lclust[-1]
        else:
            c += 1
    return clusters


def min_diff_clustering(x,
                        min_x_diff = 1.1):
    return np.fromiter(max_x_diff_iterator(x, min_x_diff),
                       dtype=int,
                       count=len(x))
