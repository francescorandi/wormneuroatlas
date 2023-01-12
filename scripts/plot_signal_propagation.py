import matplotlib.pyplot as plt, numpy as np
import wormneuroatlas as wa

# Create an instance of the NeuroAtlas
watlas = wa.NeuroAtlas()

# Get the delta F/F of the responses 
dff = watlas.get_signal_propagation_map(strain="wt")
print(np.any(np.isnan(dff)))
# or, as a shorthand alias 
# dff = watlas.get_sigpropmap(strain="wt")

# Also get the q values. 
q = watlas.get_signal_propagation_q(strain="wt")

# Reduce to head neurons.
dff = watlas.reduce_to_head(dff)
q = watlas.reduce_to_head(q)

# Plot the signal propagation map, and use the q values as transparency.
# To use the q values as transparency, substitute the nans with 1 (worst q 
# value, which will be rendered as transparent).
q[np.isnan(q)]=1.
# Make the plot. The function below does a plt.imshow(), adds neuron ids as tick
# labels, and deals with the transparency. As transparency (alphas) use 1-q,
# so that low q will be less transparent.
watlas.plot_matrix(dff,ids=watlas.head_ids,alphas=(1-q),
                   alphamin=0.4,alphamax=1.0)

plt.show()
