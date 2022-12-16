import numpy as np, h5py

data_folder = "/home/francesco/dev/wormneuroatlas/wormneuroatlas/data/"
fnames = ["cengen_021821_liberal_threshold1.csv",
          "cengen_021821_medium_threshold2.csv",
          "cengen_021821_conservative_threshold3.csv",
          "cengen_021821_stringent_threshold4.csv"]

dst = "/home/francesco/dev/wormneuroatlas/wormneuroatlas/data/cengen.h5"
f = h5py.File(dst,"w")

for i in np.arange(len(fnames)):
    data = np.genfromtxt(data_folder+fnames[i],dtype=str,delimiter=",")
    
    if i==0:
        neuron_ids = np.char.replace(data[0,3:],'"','').tolist()
        f["neuron_ids"] = neuron_ids
        
    gene_names = np.char.replace(data[1:,1],'"','').tolist()
    gene_wbids = np.char.replace(data[1:,2],'"','').tolist()
    f["gene_names_th"+str(i+1)] = gene_names
    f["gene_wbids_th"+str(i+1)] = gene_wbids

    exp = np.copy(data[1:,3:].astype(np.float32).T)
    f["tpm_th"+str(i+1)] = exp

f.close()
