import h5py, numpy as np, matplotlib.pyplot as plt

fname = "/home/francesco/dev/wormneuroatlas/wormneuroatlas/data/cengen.h5"
f = h5py.File(fname,"r+")

nids = np.array([a.decode("utf-8") for a in f["neuron_ids"][:]])
nids2 = np.array([a.decode("utf-8") for a in f["cengen_sc_1_neuron_ids"][:]])

nmatch = -1*np.ones(len(nids),dtype=int)
nfound = np.zeros(len(nids2),dtype=bool)
print("What is happening with VD_DD in the TPM and then VD and DD separately "+\
      "the thresholded matrices??") #FIXME
for ni in np.arange(len(nids)):
    m = np.where(nids[ni]==nids2)[0]
    if nids[ni] in ["AWC_OFF","AWC_ON"]: m = np.where(nids2=="AWC")[0]
    if nids[ni] in ["IL2_DV","IL2_LR"]: m = np.where(nids2=="IL2")[0]
    if nids[ni] in ["RMD_DV","RMD_LR"]: m = np.where(nids2=="RMD")[0]
    if nids[ni] in ["RME_DV","RME_LR"]: m = np.where(nids2=="RME")[0]
    if nids[ni] == "VD_DD": 
        m = np.where(nids2=="DD")[0]
        m2 = np.where(nids2=="VD")[0]
        nfound[m2] = True
    if len(m)>0:
        nmatch[ni] = m
        nfound[m] = True
print(nids2[~nfound])
        
gids = np.array([a.decode("utf-8") for a in f["gene_wbids"][:]])
gids2 = np.array([a.decode("utf-8") for a in f["cengen_sc_1_gene_wbids"][:]])

gmatch = -1*np.ones(len(gids),dtype=int)
for gi in np.arange(len(gids)):
    k = np.where(gids[gi]==gids2)[0]
    if len(k)>0:
        gmatch[gi] = k
        
cengen_sc_1b = (f["cengen_sc_1"][:][nmatch][:,gmatch]).astype(bool)
cengen_sc_2b = (f["cengen_sc_2"][:][nmatch][:,gmatch]).astype(bool)
cengen_sc_3b = (f["cengen_sc_3"][:][nmatch][:,gmatch]).astype(bool)
cengen_sc_4b = (f["cengen_sc_4"][:][nmatch][:,gmatch]).astype(bool)

##################
f["cengen_sc_1b"] = cengen_sc_1b
f["cengen_sc_2b"] = cengen_sc_2b
f["cengen_sc_3b"] = cengen_sc_3b
f["cengen_sc_4b"] = cengen_sc_4b

f.close()
