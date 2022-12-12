# wormneuroatlas
Neural signal propagation atlas (Randi et al.), genome (WormBase), and single-cell transcriptome (Taylor et al.) all in one place.

## NeuroAtlas
NeuroAtlas is the main class that aggregates all the datasets managing the conversions between different conventions for neural IDs. NeuroAtlas can be instantiated to maintain the "exact" neural identities, or to merge neurons into classes (i.e. to approximate neuron identities): `merge_bilateral=True` will merge results for, e.g., AVAL and AVAR into the class AVA_, `merge_dorsoventral=True` will merge RMED and RMEV into RME_, while `merge_numbered=True` will merge VB3, VB4, ... into VB. These options can be combined to merge, for example, SMBVL, SMBVR, SMBDL, and SMBDR into SMB__ with `merge_bilateral=True, merge_dorsoventral=True`. 
Important methods regarding neuron identities are

* NeuroAtlas.ids_to_ai() converts neuron identities to atlas-indices (ai), given the neuron-class-merging options with which the object has been created.

NeuroAtlas will use the other classes of the Python module (WormFunctionalConnectivity, WormBase, Cengen, ...) to access data from other datasets like neural signal propagation, genome, single-cell transcriptome and aggregate those datasets. 

(More details to come)
