import numpy as np, matplotlib.pyplot as plt
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
np_esconn = np.zeros((watlas.neuron_n,watlas.neuron_n))
for j in np.arange(watlas.neuron_n):
    # Find the peptides expressed in the upstream neuron j
    pep_expressed = peptide_names[pep_exp[j]!=0]
    # Find the GPCRs that bind to those peptides
    gpcrs_ = watlas.pepgpcr.get_gpcrs_binding_to(pep_expressed)
    gpcrs = [g for g_ in gpcrs_ for g in g_]
    
    
    for i in np.arange(watlas.neuron_n):
        # Find the GPCR expressed in the downstream neuron i
        gpcr_expressed = gpcr_names[gpcr_exp[i]!=0]
        
        # Count the number of GPCRs expressed in the downstream neuron that
        # are also receptors for the peptides expressed in the upstream neuron.
        # gpcrs requires two nested iterations because it 
        # gpcrs[neuropeptide] is a list of all gpcrs binding to that 
        # neuropeptide.
        for gp in gpcrs:
            if gp in gpcr_expressed: 
                np_esconn[i,j]+=1
            
plt.imshow(np_esconn, aspect="auto")
plt.show()
