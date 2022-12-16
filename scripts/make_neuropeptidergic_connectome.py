import numpy as np, matplotlib.pyplot as plt, json
import wormneuroatlas as wa

# Create NeuroAtlas object
watlas = wa.NeuroAtlas()

# Get the list of all peptide and GPCR names contained in the deorphanization
# screen.
peptide_names = np.array(watlas.pepgpcr.get_peptide_names(trim_isoforms=True))
gpcr_names = np.array(watlas.pepgpcr.get_gpcr_names(trim_isoforms=True))

# Get the neuron-resolved expression of those peptides and GPCR
pep_exp = watlas.get_gene_expression(gene_names=peptide_names)
gpcr_exp = watlas.get_gene_expression(gene_names=gpcr_names)

# Make the extrasynaptic connectome (esconn) map for neuropeptide/GPCR
pep_esconn = np.zeros((watlas.neuron_n,watlas.neuron_n))
pep_gpcr_combo = np.empty((watlas.neuron_n,watlas.neuron_n),dtype=object)
for j in np.arange(watlas.neuron_n):
    print(str(j)+"\r",end="")
    # Find the peptides expressed in the upstream neuron j
    pep_expressed = peptide_names[pep_exp[j]!=0]
    # Find the GPCRs that bind to those peptides
    gpcrs_binding = watlas.pepgpcr.get_gpcrs_binding_to(pep_expressed)
    
    for i in np.arange(watlas.neuron_n):
        # Find the GPCR expressed in the downstream neuron i
        gpcr_expressed = gpcr_names[gpcr_exp[i]!=0]
        if watlas.neuron_ids[i]=="ADAR" and watlas.neuron_ids[j] == "RMDDR":
            print("RMDDR")
            print(pep_expressed)
            print("ADAR")
            print(gpcr_expressed)
            quit()
        
        # Count the number of GPCRs expressed in the downstream neuron that
        # are also receptors for the peptides expressed in the upstream neuron.
        # gpcrs requires two nested iterations because it 
        # gpcrs[neuropeptide] is a list of all gpcrs binding to that 
        # neuropeptide.
        pep_gpcr_combo[i,j] = []
        for pep_i in np.arange(len(gpcrs_binding)):
            for gp_i in np.arange(len(gpcrs_binding[pep_i])):
                gp = gpcrs_binding[pep_i][gp_i]
                if gp in gpcr_expressed: 
                    pep_esconn[i,j]+=1
                    pep_gpcr_combo[i,j].append(pep_expressed[pep_i]+"\t"+gp)
                    
# Save pep/gpcr combos
np.save("pep_connectome_pep_gpcr_combo.npy",pep_gpcr_combo,allow_pickle=True)
print(pep_gpcr_combo[41,27])

plt.imshow(pep_esconn, aspect="auto")
plt.show()
