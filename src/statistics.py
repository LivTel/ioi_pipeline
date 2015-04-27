import numpy as np

def tukey_biweight(data, axis=None):
    data = np.array(data, copy=False)
    median = np.median(data, axis=axis)
    if axis is not None:
        delta = data - np.expand_dims(median, axis=axis)
    else:
        delta = data - median
    d = np.median(np.abs(delta), axis=axis)
    if axis is not None:
        indices = d == 0
        d[indices] = 1  # dummy value
        if np.all(indices):
            return median
        if axis is not None:
            d = np.expand_dims(d, axis=axis)
        weight = 1 - (delta / (6*d))**2
    else:
        if d == 0:
            return median
        weight = 1 - (delta / (6*d))**2
    weight = np.maximum(0, weight)
    weight = weight * weight
    if axis is not None:
        c = np.empty(indices.shape)
        c[indices] = median[indices]
        c[~indices] = (median + np.nansum(weight * delta) /
                       np.nansum(weight))[~indices]
    else:
        c = median + np.nansum(weight * delta) / np.nansum(weight)
    return c