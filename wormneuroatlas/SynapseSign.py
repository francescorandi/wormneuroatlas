import numpy as np, pandas as pd
import wormneuroatlas as wa

class SynapseSign:
    
    fname = "journal.pcbi.1007974.s003.xlsx"
    module_folder = ""
    '''Folder of the wormneuroatlas module'''
    version =  "https://doi.org/10.1371/journal.pcbi.1007974.s003|"
    description = "Container for predictions of synapse sign. "+\
                  "Currently using results from Fenyves et al. 2020."
    
    def __init__(self):
        if "\\" in wa.__file__:
            self.folder_sep = char = "\\"
        else:
            self.folder_sep = char = "/"
        self.module_folder = char.join(wa.__file__.split(char)[:-1])+char+\
                             "data"+char
                             
        self.load_from_file()
        
    def get_metadata(self):
        d = {"version": self.fname,
             "description": self.description,
             }
        
        return d
    
    def load_from_file(self):
        # Neurotransmitters first
        nt_xls = np.array(pd.read_excel(self.module_folder+self.fname, 
                                        sheet_name='1. NT expr'),
                          dtype=np.str_)
        
        self.neuron_ids = nt_xls.T[0]
        # Neuron IDs like VA01 need to be modified into VA1 to be compatible
        # with the NeuroAtlas. AWCR and AWCL need to be converted to the same
        # convention as for the connectome. Don't print a message about this
        # because these IDs only affect neurotransmitter production, and 
        # both AWCs produce Glu, so it has no effect.
        
        for i in np.arange(len(self.neuron_ids)):
            if self.neuron_ids[i][-2:].isnumeric():
                if int(self.neuron_ids[i][-2:])<10:
                    self.neuron_ids[i] = self.neuron_ids[i][:-2]+\
                                         str(int(self.neuron_ids[i][-2:]))
                                         
            if self.neuron_ids[i] == "AWCL": self.neuron_ids[i] = "AWCOFF"
            if self.neuron_ids[i] == "AWCR": self.neuron_ids[i] = "AWCON"
        
        self.dominant_nt = nt_xls.T[1]
        self.alternative_nt = nt_xls.T[2]
        
        # Receptors second
        r_xls = np.array(pd.read_excel(self.module_folder+self.fname, 
                                       sheet_name='2. Receptor gene table'),
                         dtype=np.str_)
                         
        self.neurotransmitters = ["Glu","ACh","GABA"]                 
        
        self.receptors = {}
        self.receptors["Glu"] = {}
        self.receptors["ACh"] = {}
        self.receptors["GABA"] = {}
        
        self.receptors["Glu"]["+"] = np.delete(r_xls.T[0],r_xls.T[0]=="nan")
        self.receptors["Glu"]["-"] = np.delete(r_xls.T[1],r_xls.T[1]=="nan")
        self.receptors["ACh"]["+"] = np.delete(r_xls.T[2],r_xls.T[2]=="nan")
        self.receptors["ACh"]["-"] = np.delete(r_xls.T[3],r_xls.T[3]=="nan")
        self.receptors["GABA"]["+"] = np.delete(r_xls.T[4],r_xls.T[4]=="nan")
        self.receptors["GABA"]["-"] = np.delete(r_xls.T[5],r_xls.T[5]=="nan")
        

    def get_neurons_producing(self,nt,mode="dominant"):
        '''Returns the neurons producing a given neurotransmitter.
        
        Parameters
        ----------
        nt: str
            Neurotransmitter. Can be Glu, ACh, or GABA.
        mode: str (optional)
            If dominant, the function returns the neurons for which nt is the
            dominant neurotransmitter. If alternative, only the neurons for 
            which nt is an alternative neurotransmitter. If both, all the 
            neurons producing nt in some capacity. Defalt: dominant.
        
        Returns
        -------
        neu_ids: numpy.ndarray of numpy.str_
            Neurons producing nt. Note that the ordering might be different from
            the ordering in other datasets.
        '''
        
        if mode=="dominant":
            neu_ids = self.neuron_ids[self.dominant_nt==nt]
        elif mode=="alternative":
            neu_ids = self.neuron_ids[self.alternative_nt==nt]
        elif mode=="both":
            neu_ids = self.neuron_ids[np.logical_or(self.dominant_nt==nt,
                                                    self.alternative_nt==nt)]
        
        return neu_ids
        
    def get_receptors_for(self,nt,sign):
        '''Returns the list of excitatory or inhibitory receptors for a given 
        neurotransmitter.
        
        Parameters
        ----------
        nt: str
            Neurotransmitter for which to look for receptors. Can be any of
            SynapseSign.get_neurotransmitters().
        sign: str or int
            Whether to return excitatory or inhibitory receptors. pos, exc,
            excitatory, and +1 all select excitatory receptors. neg, inh,
            inhibitory, and -1 select inhibitory receptors.
            
        Returns
        -------
        receptors: numpy.ndarray of str
            Gene names of the receptors. 
        '''
        
        if nt not in self.get_neurotransmitters():
            raise ValueError("Invalid neurotransmitter. Only the following "+\
                             "neurotransmitters are available: "+\
                             self.get_neurotransmitters())
        
        if sign in ["pos","exc","excitatory",+1]:
            sign="+"
        if sign in ["neg","inh","inhibitory",-1]: 
            sign="-"
            
        receptors = self.receptors[nt][sign]
        
        return receptors
        
    def get_neurotransmitters(self):
        return self.neurotransmitters
