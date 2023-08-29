from tools.nc_reader import nc_reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


mpl.rcParams['text.usetex'] = True

mount_loc = '/home/matthias/Documents/mount/archer2/epic3d-paper-runs/rayleigh-taylor/'


files_ref = [
    'coarsened/crse_32_epic_rt_384_fields.nc',
    'coarsened/crse_48_epic_rt_384_fields.nc',
    'coarsened/crse_64_epic_rt_384_fields.nc',
    'coarsened/crse_96_epic_rt_384_fields.nc',
    'coarsened/crse_128_epic_rt_384_fields.nc'
]

files = [
    'epic_rt_32_vmin_40_fields.nc',
    'epic_rt_48_vmin_40_fields.nc',
    'epic_rt_64_vmin_40_fields.nc',
    'epic_rt_96_vmin_40_fields.nc',
    'epic_rt_128_vmin_40_fields.nc'
    ]

labels = [
    r'$32^3$',
    r'$48^3$'
    r'$64^3$',
    r'$96^3$',
    r'$128^3$'
]

name = 'buoyancy'

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

def get_rms(x):
    return np.sqrt(np.mean(x ** 2))

def get_abse(x):
    return abs(x).max()

nc_ref = nc_reader()
nc_run = nc_reader()



try:
    for k, f in enumerate(files):
        nc_ref.open(os.path.join(mount_loc, files_ref[k]))
        tref = nc_ref.get_all('t')

        print ("Processing:", f)
        nc_run.open(os.path.join(mount_loc, f))

        ts = nc_run.get_all('t')
        n_steps = nc_run.get_num_steps() - 1

        rms = np.zeros(n_steps)
        abse = np.zeros(n_steps)

        for n in range(n_steps):

            print(n, n_steps)

            # get tref[i] <= t[step] <= tref[j]
            i, j = find_bounds(ts[n], tref)

            # we do not take the periodic edges twice into account (copy_periodic=False)
            dset_i = nc_ref.get_dataset(step=i, name=name, copy_periodic=False)
            dset_j = nc_ref.get_dataset(step=j, name=name, copy_periodic=False)

            dset_ref = time_interpl(ts[n], dset_i, tref[i], dset_j, tref[j])

            dset = nc_run.get_dataset(step=n, name=name, copy_periodic=False)

            rms[n] = get_rms(dset_ref - dset)
            abse[n] = get_abse(dset_ref - dset)

            #print(n, abse[n])

        nc_run.close()
        nc_ref.close()


        plt.plot(ts[0:n_steps], rms, label=labels[k])



    plt.xlabel('time')
    plt.ylabel('abs. max. buoyancy error')
#    plt.ylabel('relative r.m.s. buoyancy error')
    plt.savefig('abse.png', bbox_inches='tight')

except Exception as err:
    print(err)

