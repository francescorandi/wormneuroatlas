import numpy as np, matplotlib.pyplot as plt, json
import wormneuroatlas as wa

# Create NeuroAtlas object
watlas = wa.NeuroAtlas()

# This script will generate the neuropeptidergic extrasynaptic connectome 
# (pep_esconn), which can be combined with the monoaminergic connectome from
# Bentley et al. 2016, which is
ma_esconn = watlas.get_monoaminergic_connectome()

# Supplement the peptide/GPCR dataset with previous literature. Not done 
# automatically because older literature is not as systematic in measuring the
# EC50.
watlas.pepgpcr.supplement()

# Get the list of all peptide and GPCR names contained in the deorphanization
# screen.
peptide_names = np.array(watlas.pepgpcr.get_peptide_names(trim_isoforms=True))
gpcr_names = np.array(watlas.pepgpcr.get_gpcr_names(trim_isoforms=True))
gpcr_seq_ids = np.array(watlas.pepgpcr.get_gpcr_seq_id_from_name(gpcr_names))

# Get the neuron-resolved expression of those peptides and GPCR. For the GPCRs,
# supplement the names with the sequence IDs, as some entries use the sequence
# ID instead of the name.
pep_exp = watlas.get_gene_expression(gene_names=peptide_names)
gpcr_exp = watlas.get_gene_expression(gene_names=gpcr_names,
                                      gene_seq_ids=gpcr_seq_ids)
                                      
# Make the extrasynaptic connectome (esconn) map for neuropeptide/GPCR
pep_esconn = np.zeros((watlas.neuron_n,watlas.neuron_n))
# Store also the peptide/GPCR combinations (use pep_connectome_extract_combos.py
# to extract the combinations for selected pairs of neurons).
pep_gpcr_combo = np.empty((watlas.neuron_n,watlas.neuron_n),dtype=object)

for j in np.arange(watlas.neuron_n):
    # Find the peptides expressed in the upstream neuron j
    pep_expressed = peptide_names[pep_exp[j]!=0]
    # Find the GPCRs that bind to those peptides. Get the sequence IDs instead
    # of the names.
    gpcr_binding_seq_ids = \
          watlas.pepgpcr.get_gpcrs_binding_to(pep_expressed,return_seq_ids=True)
    
    for i in np.arange(watlas.neuron_n):
        # Find the GPCRs expressed in the downstream neuron i
        gpcr_expressed = gpcr_seq_ids[gpcr_exp[i]!=0]
        
        # Count the number of GPCRs expressed in the downstream neuron that
        # are also receptors for the peptides expressed in the upstream neuron.
        # gpcrs requires two nested iterations because it 
        # gpcrs[neuropeptide] is a list of all gpcrs binding to that 
        # neuropeptide.
        pep_gpcr_combo[i,j] = []
        for pep_i in np.arange(len(gpcr_binding_seq_ids)):
            for gp_i in np.arange(len(gpcr_binding_seq_ids[pep_i])):
                gp = gpcr_binding_seq_ids[pep_i][gp_i]
                if gp in gpcr_expressed: 
                    pep_esconn[i,j]+=1
                    pep_gpcr_combo[i,j].append(pep_expressed[pep_i]+"\t"+gp)
        
        # Visualize how many neurons have been processed.
        print(str(j)+"\r",end="")
                    
# Save pep/gpcr combos
np.save("pep_connectome_pep_gpcr_combo.npy",pep_gpcr_combo,allow_pickle=True)

plt.imshow(pep_esconn, aspect="auto")
plt.show()
