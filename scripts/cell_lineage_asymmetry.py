import wormneuroatlas
import matplotlib.pyplot as plt, numpy as np

atlas = wormneuroatlas.NeuroAtlas()
dff = atlas.get_signal_propagation_map()
cl = wormneuroatlas.CellLineage()

comp = atlas.get_bilateral_companions()
lin_sim = np.zeros(len(comp))*np.nan
lin_sim_measured_pairs = np.zeros(len(comp))*np.nan
done = [] # to avoid double countings

for i in range(len(comp)):
    if i not in done and comp[i] != -1:
        cell0 = atlas.neuron_ids[i]
        cell1 = atlas.neuron_ids[comp[i]]

        ls = cl.lineage_similarity(cell0,cell1)
        if ls is not None: 
            lin_sim[i] = ls
            #print(cell0,cell1,round(ls,2))
            
            if not np.isnan(dff[i,comp[i]]): lin_sim_measured_pairs[i] = ls

# Use ASEL/R as a reference, as it is known that the asymmetry in the ASEL/R
# pair is related to the 
ASEL = atlas.ids_to_ais("ASEL")

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.hist(lin_sim,bins=10,label="all bilateral",alpha=0.5,density=True)
ax.hist(lin_sim_measured_pairs,bins=10,label="all bilateral (measured)",alpha=0.5,density=True)
ax.axvline(lin_sim[ASEL],c="k",label="ASEL/R")

ax.set_xlabel("relative lineage similarity from end")
ax.set_ylabel("density")
ax.legend()

plt.show()
