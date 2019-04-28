import numpy as np
def size_factors(counts):
    counts = counts[np.alltrue(counts, axis=1)]
    logcounts = np.log(counts)
    loggeommeans = np.mean(logcounts, axis=1).reshape(len(logcounts), 1)
    sf = np.exp(np.median(logcounts - loggeommeans, axis=0))
    return sf
