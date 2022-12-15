import numpy as np, wormneuroatlas as wa

class PeptideGPCR:
    module_folder = "/".join(wa.__file__.split("/")[:-1])+"/data/"
    '''Folder of the pumpprobe module'''
    deorphan_fname = "deorphanization_media_6.csv"
    
    def __init__(self):
        data = np.loadtxt(self.module_folder+self.deorphan_fname, dtype=object,
                          skiprows=1,delimiter=",",usecols=(0,1,2,3,4,5))
        self.gpcr_seq_ids = data[:,0].astype(str)
        self.gpcr_names = data[:,1].astype(str)
        self.peptide_names = np.char.replace(data[:,2].astype(str),"pyro","")
        self.log_ec50 = data[:,3].astype(float) # in moles
        
    def get_peptide_names(self,trim_isoforms=True):
        '''Get all peptide names, potentially after trimming isoform labeling.
        
        Parameters
        ----------
        trim_isoforms: bool (optional)
            Whether to trim the isoform label from the peptide name. 
            Default: True.
        '''
        if not trim_isoforms:
            return self.peptide_names
        else:
            pep_names = self.trim_isoforms(self.peptide_names)
            return pep_names
            
    def get_gpcr_names(self,trim_isoforms=True):
        '''Get all peptide names, potentially after trimming isoform labeling.
        
        Parameters
        ----------
        trim_isoforms: bool (optional)
            Whether to trim the isoform label from the gpcr name. 
            Default: True.
        '''
        if not trim_isoforms:
            return self.gpcr_names
        else:
            gpcr_names = self.trim_isoforms(self.gpcr_names)
            return gpcr_names
        
    def get_gpcrs_binding_to(self,peptides,
                             return_seq_ids=False,trim_isoforms=True):
        '''Get the GPCRs that bind to given peptides. 
        
        Parameters
        ----------
        peptides: str or array_like of str
            Peptides for which to find GPCRs. Can be an individual peptide or a
            list of peptides. 
        return_seq_ids: bool (optional)
            Whether to return the WormBase sequence IDs instead of the protein
            names. Default: False.
        trim_isoforms: bool (optional)
            Whether to trim isoform labels from the GPCRs. Default: True.
            
        Returns
        -------
        gpcrs: str or array_like of str
            GPCRs binding to peptides. Whether a str or array_like of str 
            depends on format of input peptides.
        '''
        if type(peptides)==str:
            peptides = [peptides]
            was_scalar = True
        else:
            was_scalar = False
            
        if not return_seq_ids:
            out_info = self.gpcr_names
        else:
            out_info = self.gpcr_seq_ids
        
        gpcrs = np.empty(len(peptides),dtype=object)
        for pep_i in np.arange(len(peptides)):
            pep = peptides[pep_i].upper()
            # Determine if there is just, e.g., FLP-1 or FLP-1-1 etc.
            # In the latter case, add a dash at the end of pep, otherwise the 
            # function below will return matches also for FLP-11 etc.
            pep_ = pep+"-"
            matches_ = np.flatnonzero(np.char.find(self.peptide_names,pep_)!=-1)
                                            
            if len(matches_)>0:
                # There are FLP-1-1, etc. 
                matches = matches_
            else:
                # There is only FLP-1. Look for that exact match
                matches = self.peptide_names == pep
            
            gpcrs[pep_i] = out_info[matches]
            
            if trim_isoforms:
                # Probably faster:
                # Use regular expression to cut out the last "-1"
                gpcrs[pep_i] = self.trim_isoforms(gpcrs[pep_i])
                #gpcrs[pep_i]=["-".join(a.split("-")[:-1]) for a in gpcrs[pep_i]]
            
            gpcrs[pep_i] = np.unique(gpcrs[pep_i])
            
        if was_scalar: gpcrs = gpcrs[0]
        return gpcrs
        
    def get_peptides_binding_to(self,gpcrs,trim_isoforms=True):
        '''Get the peptides that bind to given GPCRs. 
        
        Parameters
        ----------
        gpcrs: str or array_like of str
            GPCRs for which to find peptides. Can be an individual GPCR or a
            list of GPCRs. 
        trim_isoforms: bool (optional)
            Whether to trim isoform labels from the peptides. Default: True.
            
        Returns
        -------
        peptides: str or array_like of str
            Peptides binding to GPCRs. Whether a str or array_like of str 
            depends on format of input gpcrs.
        '''
        if type(gpcrs)==str:
            gpcrs = [gpcrs]
            was_scalar = True
        else:
            was_scalar = False
            
        peptides = np.empty(len(gpcrs),dtype=object)
        for gp_i in np.arange(len(gpcrs)):
            gp = gpcrs[gp_i].upper()
            # Determine if there is just, e.g., NPR-1 or NPR-1-1 etc.
            # In the latter case, add a dash at the end of pep, otherwise the 
            # function below will return matches also for NPR-11 etc.
            gp_ = gp+"-"
            matches_ = np.flatnonzero(np.char.find(self.gpcr_names,gp_)!=-1)
            
            if len(matches_)>0:
                # There are FLP-1-1, etc. 
                matches = matches_
            else:
                # There is only FLP-1. Look for that exact match
                matches = self.gpcr_names == gp
                
            peptides[gp_i] = self.peptide_names[matches]
            
            if trim_isoforms:
                # Probably faster:
                # Use regular expression to cut out the last "-1"
                peptides[gp_i] = self.trim_isoforms(peptides[gp_i])
                #peptides[gp_i]=["-".join(a.split("-")[:-1]) \
                #                    for a in peptides[gp_i]]
            
            peptides[gp_i] = np.unique(peptides[gp_i])
            
            if was_scalar: peptides = peptides[0]
            return peptides
   
    @staticmethod     
    def trim_isoforms(b,n=2):
        '''Trim isoform label ("-n") from a protein name.
        
        Parameters
        ----------
        b: str
            Protein name to be trimmed.
        n: int (optional)
            How many dash-separated parts of the string to expect in a 
            proteint name without any isoform labeling (e.g. n=2 for FLP-1).
        
        Returns
        -------
        b_: str
            Trimmed protein name.
        '''
        if type(b)==str:
            b = [b]
            was_scalar = True
        else:
            was_scalar = False
        
        b_ = []
        for a in b:
            if len(a.split("-"))>n:
                b_.append("-".join(a.split("-")[:-1]))
        
        if was_scalar: b_ = b_[0]
        return b_
