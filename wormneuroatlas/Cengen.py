import numpy as np, h5py,matplotlib.pyplot as plt
import wormneuroatlas as wa

class Cengen:
    
    module_folder = ""
    '''Folder of the wormneuroatlas module'''
    cengen_fname = "cengen.h5"
    
    def __init__(self,wormbase=None):
        if "\\" in wa.__file__:
            self.folder_sep = char = "\\"
        else:
            self.folder_sep = char = "/"
        self.module_folder = char.join(wa.__file__.split(char)[:-1])+char+\
                             "data"+char
        
        self.h5 = h5py.File(self.module_folder+self.cengen_fname,"r")
        # Directly load neuron from file. These are the same for each 
        # threshold file.
        self.neuron_ids = [n.decode("utf-8") for n in self.h5["neuron_ids"][:]]
        self.neuron_ids = np.array(self.neuron_ids)
        self.neuron_n = self.n_neurons = len(self.neuron_ids)
        
        # These are variable that come with different gene orderings in the 
        # external datasets.
        self.gene_expression = np.zeros(5,dtype=object)
        self.gene_wbids = np.zeros(5,dtype=object)
        self.gene_wbids_int = np.zeros(5,dtype=object)
        self.gene_names = np.zeros(5,dtype=object)
        self.gene_n = self.n_genes = np.zeros(5,dtype=int)
        
        for th in np.arange(5)[1:]:
            # Load gene IDs from file.
            a=[g.decode("utf-8") for g in self.h5["gene_wbids_th"+str(th)][:]]
            self.gene_wbids[th] = np.array(a)
            # Make a integer version of the WormBase gene IDs.
            self.gene_wbids_int=np.array([int(a[6:]) for a in a])
            # Load gene names from file.
            a = [g.decode("utf-8") for g in self.h5["gene_names_th"+str(th)][:]]
            self.gene_names[th] = np.array(a)
            # Store the number of genes
            self.gene_n[th] = len(a)
            
            self.gene_expression[th] = self.h5["tpm_th"+str(th)][:]
            
        # Store the instance of wormbase. Only necessary for specific things,
        # like looking up gene_wbids from gene_names.
        self.wormbase = wormbase
        
    def get_neuron_ids(self):
        return self.neuron_ids
        
    def get_gene_wbids(self,th):
        return self.gene_wbids[th]
        
    def get_gene_names(self,th):
        return self.gene_names[th]
        
    def get_expression(self, gene_cis=None, gene_wbids=None, 
                       gene_names=None, gene_seq_ids=None,
                       neuron_cis=None, neuron_ids=None,
                       th=4):
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
        gene_seq_ids: array_like of strings (optional)
            Gene sequence IDs. Only used if gene_names is not None. If a gene
            name is unknown, an attempt is made to match the sequence ID. 
            Default: None.
        neuron_cis: "all" or array_like of int (optional)
            CeNGEN indices of the neurons. If None or "all", all are selected.
            Default: None.
        neuron_ids: array_like of strings (optional)
            CeNGEN-style neuron IDs. Only used if neuron_cis is None. Default:
            None.
        th: int (optional)
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
        gene_cis = self.parse_genes_input(th,gene_cis,gene_wbids,gene_names)
        unknown_genes = gene_cis==-1
        gene_cis[gene_cis==-1] = 0

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
        else:
            neuron_cis = np.arange(self.neuron_n)
                
        neuron_cis = np.array(neuron_cis)
        
        expression = self.gene_expression[th][neuron_cis][:,gene_cis]
        
        # Reapply mask for unknown genes
        expression[:,unknown_genes] = 0
        
        return expression
        
    
    def parse_genes_input(self,th=4,
                          gene_cis=None, gene_wbids=None, 
                          gene_names=None, gene_seq_ids=None):
        '''Parses flexible input of gene selection.
        
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
        gene_seq_ids: array_like of strings (optional)
            Gene sequence IDs. If gene_names is not None and a gene name is 
            unknown, an attempt is made to match the sequence ID. This is 
            necessary because some genes are labeled with their sequence ID
            and not with a name. If gene_cis, gene_wbids, and gene_names are 
            None, then gene_seq_ids are used independently to query for genes. 
            Default: None.
            
        Returns
        -------
        gene_cis: numpy.ndarray of int
            CeNGEN indices of the genes.
        '''
        if gene_cis is not None:
            # Allow for "all" or a specific set
            if type(gene_cis)==str:
                if gene_cis=="all":
                    gene_cis = np.arange(self.gene_n[th])
            else:
                # You already have the CeNGEN indices of the genes
                pass
                
        elif gene_wbids is not None:
            # If the WormBase IDs are passed, they could be the integer or
            # string version of the ID (WBGene00000001 or 1)
            gene_cis = []
            for gwbid in gene_wbids:
                if type(gwbid)==int:
                    gene_cis.append(np.where(gwbid==self.gene_wbids_int[th])[0][0])
                elif type(gwbid)==str:
                    gene_cis.append(np.where(gwbid==self.gene_wbids[th])[0][0])
                else:
                    gene_cis.append(None)
                    
        elif gene_names is not None:
            gene_names = [g.lower() for g in gene_names]
            # If names are passed, use the WormBase instance to transform those
            # names to WormBase IDs and then transform the latter to CeNGEN IDs.
            gene_cis = []
            for iname in np.arange(len(gene_names)):
                gname = gene_names[iname]
                match = np.where(gname==self.gene_names[th])[0]
                if len(match)>0:
                    gene_cis.append(match[0])
                else:
                    if gene_seq_ids is not None:
                        seq_id = gene_seq_ids[iname]
                        match = np.where(seq_id==self.gene_names[th])[0]
                        if len(match)>0:
                            gene_cis.append(match[0])
                        else:
                            gene_cis.append(-1)
                    else:
                        gene_cis.append(-1)
                        
        elif gene_seq_ids is not None:
            gene_cis = []
            for iseq in np.arange(len(gene_seq_ids)):
                seq_id = gene_seq_ids[iseq]
                match = np.where(seq_id==self.gene_seq_ids[th])[0]
                if len(match)>0:
                    gene_cis.append(match[0])
                else:
                    gene_cis.append(-1)

        else:
            #If nothing is specified, assume all genes are selected.
            gene_cis = np.arange(self.gene_n[th])
        
        gene_cis = np.array(gene_cis)   
                    
        return gene_cis
        
    
    @staticmethod
    def slice_h5_2d(h5,slice1,slice2):
        '''Slice h5 dealing with the fact that the indices need to be sorted.
        Sort and de-sort.
        '''

        # This is for multiple 0s.
        slice1,indices1 = np.unique(slice1,return_inverse=True)
        slice2,indices2 = np.unique(slice2,return_inverse=True)
        
        is_slice1_sorted = np.all(slice1[:-1] <= slice1[1:])
        is_slice2_sorted = np.all(slice2[:-1] <= slice2[1:])
        
        if is_slice1_sorted:
            slice1_ = slice1
            sorter1 = inv_sorter1 = np.arange(len(slice1))
        else:
            sorter1 = np.argsort(slice1) 
            inv_sorter1 = np.empty_like(sorter1)
            inv_sorter1[sorter1] = np.arange(len(sorter1))
            slice1_ = slice1[sorter1]
            
        if is_slice2_sorted:
            slice2_ = slice2
            sorter2 = inv_sorter2 = np.arange(len(slice2))
        else:
            sorter2 = np.argsort(slice2) 
            inv_sorter2 = np.empty_like(sorter2)
            inv_sorter2[sorter2] = np.arange(len(sorter2))
            slice2_ = slice2[sorter2]
            
        h5_sliced = h5[slice1_][:,slice2_]
        h5_sliced = h5_sliced[inv_sorter1][:,inv_sorter2]
        h5_sliced = h5_sliced[indices1][:,indices2]
        
        return h5_sliced
            
            
         
                    
                    
                    
