import numpy as np

def find_closest_index(val, arr):
    # 27 July 2023
    # https://stackoverflow.com/a/2566508
    return (np.abs(arr - val)).argmin()

# find i and j such that arr[i] <= val <= arr[j]:
def find_bounds(val, arr):
    i = find_closest_index(val, arr)
    j = i + 1
    if arr[i] > val:
        j = i
        i = i - 1

    if arr[i] > val or arr[j] < val:
        raise ValueError('Could not find the correct bounds.')

    return i, j

# t1 <= t <= t2
def time_interpl(t, x1, t1, x2, t2):

    dt1 = t - t1
    dt2 = t2 - t

    dt = t2 - t1

    if dt == 0.0:
        raise ValueError('Time step difference is zero.')


    w1 = 1 - dt1 / dt
    w2 = 1 - dt2 / dt

    return w1 * x1 + w2 * x2

def get_dataset(ncr, name, t, copy_periodic=False):
    ts = ncr.get_all('t')

    # get tref[i] <= t[step] <= tref[j]
    i, j = find_bounds(t, ts)

    # we do not take the periodic edges twice into account (copy_periodic=False)
    dset_i = ncr.get_dataset(step=i, name=name, copy_periodic=copy_periodic)
    dset_j = ncr.get_dataset(step=j, name=name, copy_periodic=copy_periodic)

    return time_interpl(t, dset_i, ts[i], dset_j, ts[j])
