#!/usr/bin/python3
#!/usr/bin/env python
#-*- coding: utf-8 -*-
#==============================================================================
# Copyright (c) 2018, William Gurecky
# All rights reserved.
#
# DESCRIPTION:
# Inspired from stack overflow question:
# https://stackoverflow.com/questions/52227599/interpolate-griddata-uses-only-one-core
#
# OPTIONS: --
# AUTHOR: William Gurecky
# CONTACT: william.gurecky@gmail.com
#==============================================================================
import numpy as np
from scipy import interpolate
import dask.array as da
import matplotlib.pyplot as plt
from dask import delayed

def dask_gd2(xx, yy, z_array, target_xi, target_yi, algorithm='cubic', **kwargs):
    """!
    @brief general parallel interpolation using dask and griddata
    @param xx 1d or 2d array of x locs where data is known
    @param yy 1d or 2d array of x locs where data is known
    @param z_array 1d or 2d array of x locs where data is known
    @param target_xi 2d array (or 1d grid spacing array)
    @param target_yi 2d array (or 1d grid spacing array)
    """
    n_jobs = kwargs.pop("n_jobs", 4)
    chunk_size = kwargs.get("chunk_size", int(xx.size / (n_jobs - 1)))
    if len(target_xi.shape) < 2:
        xxt, yyt = np.meshgrid(target_xi, target_yi)
    elif len(target_xi.shape) > 2:
        raise RuntimeError
    else:
        xxt, yyt = target_xi, target_yi
    assert xxt.shape == yyt.shape
    z_target = np.full(np.shape(xxt), np.nan)

    # evenly mix nans into dataset.  nans mark where data is needed
    n_splits = n_jobs * 8
    sp_xx, sp_xxt = np.array_split(xx.flatten(), n_splits), np.array_split(xxt.flatten(), n_splits)
    sp_yy, sp_yyt = np.array_split(yy.flatten(), n_splits), np.array_split(yyt.flatten(), n_splits)
    sp_zz, sp_zzt = np.array_split(z_array.flatten(), n_splits), np.array_split(z_target.flatten(), n_splits)

    all_x = np.concatenate(np.array((sp_xx, sp_xxt)).T.flatten())
    all_y = np.concatenate(np.array((sp_yy, sp_yyt)).T.flatten())
    all_z = np.concatenate(np.array((sp_zz, sp_zzt)).T.flatten())

    # make dask arrays
    import pdb; pdb.set_trace()
    dask_xx = da.from_array(all_x, chunks=chunk_size, name="dask_x")
    dask_yy = da.from_array(all_y, chunks=chunk_size, name="dask_y")
    dask_zz = da.from_array(all_z, chunks=chunk_size, name="dask_z")

    dask_valid_x1 = dask_xx[~da.isnan(dask_zz)]
    dask_valid_y1 = dask_yy[~da.isnan(dask_zz)]
    dask_valid_z1 = dask_zz[~da.isnan(dask_zz)]

    # where to interplate to
    dask_target_x = dask_xx[da.isnan(dask_zz)]
    dask_target_y = dask_yy[da.isnan(dask_zz)]

    # interpolate for missing values
    zz_grid = dask_interpolate(dask_valid_x1, dask_valid_y1, dask_valid_z1, dask_target_x, dask_target_y, algorithm=algorithm, **kwargs)
    return zz_grid.reshape(xxt.shape)


def dask_gd2_nanfill(xx, yy, z_array, algorithm='cubic', **kwargs):
    """!
    @brief 2d interpolation using dask and griddata
    @param xx np_2darray x coord array
    @param yy np_2darray y coord array
    @param z_array np_2darray response vals
    """
    n_jobs = kwargs.pop("n_jobs", 4)
    chunk_size = kwargs.get("chunk_size", int(xx.size / (n_jobs - 1)))
    # make dask arrays
    dask_xyz = da.from_array((xx, yy, z_array), chunks=(3, chunk_size, "auto"), name="dask_all")
    dask_xx = dask_xyz[0,:,:]
    dask_yy = dask_xyz[1,:,:]
    dask_zz = dask_xyz[2,:,:]

    # select only valid values
    dask_valid_x1 = dask_xx[~da.isnan(dask_zz)]
    dask_valid_y1 = dask_yy[~da.isnan(dask_zz)]
    dask_valid_z1 = dask_zz[~da.isnan(dask_zz)]

    # interpolate for missing values
    return dask_interpolate(dask_valid_x1, dask_valid_y1, dask_valid_z1, dask_xx, dask_yy, algorithm=algorithm, **kwargs)


def dask_interpolate(dask_valid_x1, dask_valid_y1, dask_valid_z1, dask_xx, dask_yy, algorithm='cubic', vis_out='dask_par.png'):
    # gd_chunked = [delayed(rbf_wrapped)(x1, y1, newarr, xx, yy) for \
    gd_chunked = [delayed(gd_wrapped)(x1.flatten(), y1.flatten(), newarr.flatten(), xx, yy, algorithm) for \
                x1, y1, newarr, xx, yy \
                in \
                zip(dask_valid_x1.to_delayed().flatten(),
                    dask_valid_y1.to_delayed().flatten(),
                    dask_valid_z1.to_delayed().flatten(),
                    dask_xx.to_delayed().flatten(),
                    dask_yy.to_delayed().flatten())]
    gd_out = delayed(da.concatenate)(gd_chunked, axis=0)
    gd_out.visualize(vis_out)
    gd1 = np.array(gd_out.compute())
    print(gd1)
    print(gd1.shape)

    # prove we have no more nans in the data
    assert ~np.isnan(np.sum(gd1))
    return gd1


def rbf_wrapped(x1, y1, newarr, xr, yr):
    print("local x.size: ", x1.size)
    rbf_interpolant = interpolate.Rbf(x1, y1, newarr, function='linear')
    return rbf_interpolant(xr, yr)


def gd_wrapped(x, y, v, xp, yp, algorithm='cubic', extrapolate=True):
    # source: https://programtalk.com/python-examples/scipy.interpolate.griddata.ravel/
    print("local x.size: ", x.size)
    if x.size == 0 or xp.size == 0:
        # nothing to do
        return np.array([[]])
    if algorithm not in ['cubic', 'linear', 'nearest']:
        raise ValueError("Invalid interpolation algorithm: " + str(algorithm))
    # known data (x, y, v)  can be either 1d or 2d arrays of same size
    # target grid: (xp, yp), xp, yp must be 2d arrays of the same shape
    grid = interpolate.griddata((x, y), v, (xp, yp),
                                method=algorithm, rescale=True)
    if extrapolate and algorithm != 'nearest' and np.any(np.isnan(grid)):
        grid = extrapolate_nans(xp, yp, grid)
    return grid


def extrapolate_nans(x, y, v):
    if x.size == 0 or y.size == 0:
        # nothing to do
        return np.array([[]])
    if np.ma.is_masked(v):
        nans = v.mask
    else:
        nans = np.isnan(v)
    notnans = np.logical_not(nans)
    v[nans] = interpolate.griddata((x[notnans], y[notnans]), v[notnans],
                                   (x[nans], y[nans]),
                                   method='nearest', rescale=True)
    return v


if __name__ == "__main__":
    """
    # create data with random missing entries
    ar_size_x, ar_size_y = 500, 600
    chunk_size = 100
    z_array = np.ones((ar_size_x, ar_size_y))
    z_array[np.random.randint(0, ar_size_x-1, 50),
          np.random.randint(0, ar_size_y-1, 50)]= np.nan

    # XY coords
    x = np.linspace(0, 3, z_array.shape[1])
    y = np.linspace(0, 3, z_array.shape[0])

    # gen sin wave for testing
    z_array = z_array * np.sin(x)
    # prove there are nans in the dataset
    assert np.isnan(np.sum(z_array))

    xx, yy = np.meshgrid(x, y)
    print("global x.size: ", xx.size)

    # perform example parallel interp
    gd1 = dask_gd2_nanfill(xx, yy, z_array, algorithm='cubic')
    plt.figure()
    plt.imshow(gd1)
    plt.savefig("dask_par_sin_t1.png")
    plt.close()
    """

    # test2
    ar_size_x, ar_size_y = 500, 600
    chunk_size = 100
    z_array = np.ones((ar_size_x, ar_size_y))

    # XY coords
    x = np.linspace(0, 3, z_array.shape[1])
    y = np.linspace(0, 3, z_array.shape[0])
    xx, yy = np.meshgrid(x, y)

    # gen sin wave for testing
    z_array = z_array * np.sin(x)

    # target grid
    x_target = np.linspace(0, 3, 200)
    y_target = np.linspace(0, 3, 100)

    gd1 = dask_gd2(xx, yy, z_array, x_target, y_target, algorithm='cubic', n_jobs=8)
    plt.figure()
    plt.imshow(gd1)
    plt.savefig("dask_par_sin_t2.png")
    plt.close()
