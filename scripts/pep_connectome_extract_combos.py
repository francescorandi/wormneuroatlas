import sys, numpy as np
import wormneuroatlas as wa

watlas = wa.NeuroAtlas()

# Allow to pass a file containing the pairs for which to save the combo of
# neuropeptide and GPCR.
extr_combo_for_fname = sys.argv[1]

# Read in the requested pairs from the file and store the combo entries for 
# that pair.
f_extr_for = open(extr_combo_for_fname,"r")
f_extr = open(extr_combo_for_fname.split(".")[0]+"_pepgpcr_combos.txt","w")
combos = np.load("pep_connectome_pep_gpcr_combo.npy",allow_pickle=True)
for l in f_extr_for.readlines():
    neu_ids = l.rstrip().split("->")[::-1]
    pair = watlas.ids_to_ais(neu_ids)
    if neu_ids[1]=="AQR":
        print(pair)
    f_extr.write(l.rstrip()+"\n")
    for combo in combos[pair[0],pair[1]]:
        f_extr.write("\t"+combo+"\n")
    
f_extr_for.close()
f_extr.close()
