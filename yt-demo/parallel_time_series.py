import numpy as np
import matplotlib.pyplot as plt
import yt
# If you want to run in parallel this *must* be the first
# thing after importing yt!
yt.enable_parallelism()

# We first create a time series object of the data files, using a
# "wildcard" syntax

ts = yt.DatasetSeries("../data/WindTunnel/windtunnel_4lev_hdf5_plt_cnt_*")

# This is a dictionary object which we'll use to store the results in
my_storage = {}

# Now we loop over the time series, calculating the average x-velocity and
# storing the result. We also make slice plots of every snapshot,
# annotating the grid lines on top.

for sto, ds in ts.piter(storage=my_storage):
    dd = ds.all_data() # This is an object giving us all the data
    # This line computes the average x-velocity weighted by the
    # density
    vx = dd.quantities.weighted_average_quantity("velocity_x", "density")

    # We now store both the filename and
    # the result in the storage.
    sto.result_id = str(ds) # Taking str() of ds gives us the filename
    sto.result = (ds.current_time, vx)

    # Make a slice plot, setting the z-lim of the density field,
    # drawing the grid lines on top, and saving it.
    slc = yt.SlicePlot(ds, "z", ["density"])
    slc.set_zlim("density", 0.8, 8.0)
    slc.set_log("density", False)
    slc.annotate_grids()
    slc.save()

# Now, let processor 0 gather the values of the time and velocity
# and save them to a plot.
if yt.is_root():
    time = np.array([v[0] for v in my_storage.values()])
    vx = np.array([v[1] for v in my_storage.values()])
    idxs = np.argsort(time)
    plt.plot(time[idxs], vx[idxs])
    plt.savefig("vx_vs_time.png")
