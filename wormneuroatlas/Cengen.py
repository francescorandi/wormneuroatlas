import numpy as np, h5py
import wormneuroatlas as wa

class Cengen:
    
    module_folder = "/".join(wa.__file__.split("/")[:-1])+"/data/"
    '''Folder of the pumpprobe module'''
    cengen_fname = "cengen.h5"
    
    def __init__(self,wormbase=None):
        self.h5 = h5py.File(self.module_folder+self.cengen_fname,"r")
        # Directly load neuron and gene IDs from file.
        self.neuron_ids = [n.decode("utf-8") for n in self.h5["neuron_ids"][:]]
        self.neuron_ids = np.array(self.neuron_ids)
        self.neuron_n = self.n_neurons = len(self.neuron_ids)
        self.gene_wbids = [g.decode("utf-8") for g in self.h5["gene_wbids"][:]]
        self.gene_wbids = np.array(self.gene_wbids)
        self.gene_n = self.n_genes = len(self.gene_wbids)
        # Make a integer version of the WormBase gene IDs.
        self.gene_wbids_int = np.array([int(a[6:]) for a in self.gene_wbids])
        
        # Don't load gene expression from file yet. Only on demand.
        self.gene_expression = self.h5["cengen_TPM"]
        
        # Make the masks for different confidences
        self.confidence_mask = np.ones((5,self.neuron_n,self.gene_n),dtype=bool)
        for cmi in np.arange(5)[1:]:
            print("Still something to fix for the confidence thresholds.")
            #self.confidence_mask[cmi] = self.h5["cengen_sc_"+str(cmi)][:]
        
        # Store the instance of wormbase. Only necessary for specific things,
        # like looking up gene_wbids from gene_names.
        self.wormbase = wormbase
        
        
    def get_neuron_ids(self):
        return self.neuron_ids
        
    def get_expression(self, gene_cis=None, gene_wbids=None, gene_names=None,
                       neuron_cis=None, neuron_ids=None,
                       threshold=4):
        '''Get gene expression for given genes and neurons.
        
        Parameters
        ----------
        gene_cis: "all" or array_like of int (optional)
            CeNGEN indices of the genes. If None or "all", all are selected.
            Default: None.
        gene_wbids: array_like of int or strings (optional)
            WormBase IDs of the genes. Can be string or integer (WB00000001 or 
            1). Only used if gene_cis is None. Default: None.
        gene_names: array_like of strings (optional)
            Gene names. Only used if both gene_cis and gene_wbids are None.
            Default: None.
        neuron_cis: "all" or array_like of int (optional)
            CeNGEN indices of the neurons. If None or "all", all are selected.
            Default: None.
        neuron_ids: array_like of strings (optional)
            CeNGEN-style neuron IDs. Only used if neuron_cis is None. Default:
            None.
        threshold: int (optional)
            Confidence threshold, as in the CeNGEN online app. Can range from
            0 (no threshold) to 4. Default: 4.
            
        Returns
        -------
        expression: 2D numpy.ndarry
            expression[neuron_i,gene_j] is the expression level of the j-th of 
            the input genes in the the i-th neuron of the input neurons.
        '''
        
        ########################################################################
        # Parse the different possible inputs to select which genes to use.
        ########################################################################
        if gene_cis is not None:
            # Allow for "all" or a specific set
            if type(gene_cis)==str:
                if gene_cis=="all":
                    gene_cis = np.arange(self.gene_n)
            else:
                # You already have the CeNGEN indices of the genes
                pass
        elif gene_wbids is not None:
            # If the WormBase IDs are passed, they could be the integer or
            # string version of the ID (WBGene00000001 or 1)
            gene_cis = []
            for gwbid in gene_wbids:
                if type(gwbid)==int:
                    gene_cis.append(np.where(gwbid==self.gene_wbids_int)[0][0])
                elif type(gwbid)==str:
                    gene_cis.append(np.where(gwbid==self.gene_wbids)[0][0])
                else:
                    gene_cis.append(None)
                    
        elif gene_names is not None:
            # If names are passed, use the WormBase instance to transform those
            # names to WormBase IDs and then transform the latter to CeNGEN IDs.
            if self.wormbase is not None:
                gene_cis = []
                for gn in gene_names:
                    # Get WormBase ID
                    wbid = self.wormbase.name_to_gene_wbid(gn,alt=True)
                    # Get CeNGEN gene_i
                    ci = np.where(wbid==self.gene_wbids_int)[0][0]
                    gene_cis.append(ci)

        else:
            #If nothing is specified, assume all genes are selected.
            gene_cis = np.arange(self.gene_n)
        
        gene_cis = np.array(gene_cis)


        #########################            
        # Parse the neuron input.
        #########################
        if neuron_cis is not None:
            if type(neuron_cis)==str:
                if neuron_cis == "all":
                    neuron_cis = np.arange(self.neuron_n)
        elif neuron_ids is not None:
            neuron_cis = []
            for nid in neuron_ids:
                nci = np.where(nid==self.neuron_ids)[0][0]
                neuron_cis.append(nci)
                
        neuron_cis = np.array(neuron_cis)
        
        expression = self.gene_expression[neuron_cis][:,gene_cis]
        # Apply mask given the threshold on the confidence
        expression *= self.confidence_mask[threshold][neuron_cis][:,gene_cis]
        
        return expression
                    
                    
                    
