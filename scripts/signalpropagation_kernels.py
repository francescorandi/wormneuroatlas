import matplotlib.pyplot as plt, numpy as np
import wormneuroatlas as wa

# Create an instance of the NeuroAtlas
watlas = wa.NeuroAtlas()

# Make a time axis on which to evaluate the kernel. NeuroAtlas returns an 
# object representing the kernel that can be evaluated on any (all-positive)
# time axis, and not the kernel already evaluated as an array.
time = np.linspace(0,20,1000)

# Get the kernel for a specific pair. If there is no kernel in the data, None
# is returned
kernel = watlas.get_kernel(i=0,j=58,strain="wt")
if kernel is not None:
    # Evaluate the kernel on the given time-axis. 
    y = kernel.eval(time)
    # or in the short-hand notation
    y = kernel(time)
    
    plt.plot(time,y)
    plt.show()
