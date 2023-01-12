import numpy as np, re, wormneuroatlas as wa

class PeptideGPCR:
    module_folder = ""
    '''Folder of the wormneuroatlas module'''
    deorphan_fname = "deorphanization_media_6.csv"
    bentley_fname = "esconnectome_neuropeptides_Bentley_2016.csv"
    froonikcx_fname = "froonikcx_peptide_gpcr.txt"
    cached_seq_id_fname = "cached_seq_ids.txt"
    
    def __init__(self):
        if "\\" in wa.__file__:
            self.folder_sep = char = "\\"
        else:
            self.folder_sep = char = "/"
        self.module_folder = char.join(wa.__file__.split(char)[:-1])+char+\
                             "data"+char
        
        # Load data from Beets et al. 2022
        data = np.loadtxt(self.module_folder+self.deorphan_fname, dtype=object,
                          skiprows=1,delimiter=",",usecols=(0,1,2,3,4,5))
        self.gpcr_seq_ids = data[:,0].astype(str)
        self.gpcr_names = data[:,1].astype(str)
        self.peptide_names = np.char.replace(data[:,2].astype(str),"pyro","")
        for pi in np.arange(len(self.peptide_names)):
            self.peptide_names[pi]=re.sub(r' \(.*\)',"",self.peptide_names[pi])
        self.log_ec50 = data[:,3].astype(float) # log10 of moles
        
        self.cached_seq_ids = np.char.lower(np.loadtxt(
                        self.module_folder+self.cached_seq_id_fname,
                        dtype=str))
        
    def supplement(self,bentley=True,froonikcx=True,verbose=False):
        '''Supplement the peptide/GPCR combinations using additional works
        from the literature. Not done by default, because the main source 
        for the deorphanization used in the __init__() reports also 
        activation thresholds and systematically considers low-threshold 
        binding between peptides and GPCR.
        
        Parameters
        ----------
        bentley: bool (optional)
            Whether to use data from Bentley et al. 2016, Multilayer connectome
            of C. elegans. Default: True.
        froonikcs: bool (optional)
            Whether to use data from Froonikcx et al. 2012. Default: True
        verbose: bool (optional)
            Whether to print out what combinations of peptides and GPCR are
            added from these supplemental datasets. Default: False.
        '''
        
        if bentley or froonikcx:
            comb_a_ = []
            for i in np.arange(len(self.peptide_names)):
                a = self.trim_isoforms([self.peptide_names[i]])[0]
                b = self.gpcr_seq_ids[i]
                comb_a_.append(a+"->"+b)
                
        if bentley:
            # Supplement data from Bentley et al. 2016
            data = np.loadtxt(self.module_folder+self.bentley_fname,dtype=str,
                              skiprows=1,delimiter=",",usecols=(2,3))
            comb_b_ = [d[0]+"->"+d[1] for d in data]
            comb_b_ = np.unique(comb_b_)
            
            # use the array as comb_b_ with 
            # trim_isoforms(self.peptide_names) -> self.gpcr_names, 
            # and then check if comb_b_[i] in comb_a_
            if verbose: print("ADDED THROUGH BENTLEY")
            for cb_ in comb_b_:
                cb_ = cb_.upper()
                if cb_ not in comb_a_:
                    if verbose: print("\t"+cb_)
                    cb = cb_.split("->")
                    self.peptide_names = np.append(self.peptide_names,cb[0])
                    gp_seq_id = self.get_gpcr_seq_id_from_name(cb[1])
                    self.gpcr_names = np.append(self.gpcr_names,cb[1])
                    self.gpcr_seq_ids = np.append(self.gpcr_seq_ids,gp_seq_id)
                    self.log_ec50 = np.append(self.log_ec50,0)
                    
        if froonikcx:      
            # Supplement data from Froonickx (see table at the bottom)
            # https://www.frontiersin.org/articles/10.3389/fendo.2012.00167/full
            
            data = np.loadtxt(self.module_folder+self.froonikcx_fname,dtype=str)
            # Use seq_ids instead of gene names
            comb_b_ = [self.trim_isoforms(d[0])+"->"+self.trim_isoforms(d[1])
                        for d in data]
            comb_b_ = np.unique(comb_b_,)
            
            if verbose: print("ADDED THROUGH FROONIKCX")
            for i in np.arange(len(comb_b_)):
                cb_ = comb_b_[i]
                cb_ = cb_.upper()
                if cb_ not in comb_a_:
                    if verbose: print("\t"+cb_)
                    cb = cb_.split("->")
                    self.peptide_names = np.append(self.peptide_names,cb[0])
                    gp_seq_id = self.get_gpcr_seq_id_from_name(cb[1])
                    self.gpcr_names = np.append(self.gpcr_names,cb[1])
                    self.gpcr_seq_ids = np.append(self.gpcr_seq_ids,gp_seq_id)
                    self.log_ec50 = np.append(self.log_ec50,0)
                    
        
    def get_peptide_names(self,trim_isoforms=True):
        '''Get all peptide names, potentially after trimming isoform labeling.
        
        Parameters
        ----------
        trim_isoforms: bool (optional)
            Whether to trim the isoform label from the peptide name. 
            Default: True.
        '''
        if not trim_isoforms:
            return np.unique(self.peptide_names)
        else:
            pep_names = np.unique(self.trim_isoforms(self.peptide_names))
            return pep_names
            
    def get_gpcr_names(self,trim_isoforms=True):
        '''Get all GPCR names, potentially after trimming isoform labeling.
        
        Parameters
        ----------
        trim_isoforms: bool (optional)
            Whether to trim the isoform label from the gpcr name. 
            Default: True.
        '''
        if not trim_isoforms:
            gpcr_names = np.unique(self.gpcr_names)
        else:
            gpcr_names = np.unique(self.trim_isoforms(self.gpcr_names))
        
        return gpcr_names
                
    def get_gpcr_seq_ids(self,trim_isoforms=True):
        '''Get all GPCR sequence IDs, potentially after trimming isoform 
        labeling.
        
        Parameters
        ----------
        trim_isoforms: bool (optional)
            Whether to trim the isoform label from the gpcr sequence id. 
            Default: True.
        '''
        if not trim_isoforms:
            gpcr_seq_ids = np.unique(self.gpcr_seq_ids)
        else:
            gpcr_seq_ids = np.unique(self.trim_isoforms(self.gpcr_seq_ids,n=1))
                            
        return gpcr_seq_ids
        
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
                # There are FLP-1-1, etc. Use the matches of FLP-1-* BUT ALSO
                # the exact matches to FLP-1, which would otherwise be 
                # discarded.
                matches = matches_
                exact_matches = np.where(self.peptide_names == pep)[0]
                matches = np.append(matches, exact_matches)
            else:
                # There is only FLP-1. Look for that exact match
                matches = self.peptide_names == pep
            
            gpcrs[pep_i] = out_info[matches]
            
            if trim_isoforms:
                if return_seq_ids: n=1
                else: n=2
                gpcrs[pep_i] = self.trim_isoforms(gpcrs[pep_i],n=n)
            
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
                # There are NPR-1-1, etc. Use the matches of NPR-1-* BUT ALSO
                # the exact matches to NPR-1, which would otherwise be 
                # discarded.
                matches = matches_
                exact_matches = np.where(self.gpcr_names == gp)[0]
                matches = np.append(matches, exact_matches)
            else:
                # There is only NPR-1. Look for that exact match
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
            
    def get_gpcr_seq_id_from_name(self,names):
        if type(names) in [str,np.str_]: 
            names=[names]
            was_scalar = True
        else:
            was_scalar = False
        
        seq_ids = []
        for name in names:
            name_ = name+"-"
            matches_ = np.flatnonzero(np.char.find(self.gpcr_names,name_)!=-1)
            
            if len(matches_)>0:
                # There are NPR-1-1, etc. Use the matches of NPR-1-* BUT ALSO
                # the exact matches to NPR-1, which would otherwise be 
                # discarded.
                matches = matches_
                
                # No matter if there are NPR-1-1, if the name you passed has no
                # isoform label, then any of the NPR-1-* will give you the
                # sequence name you want, after trimming.
                seq_id = self.gpcr_seq_ids[matches[0]]
                seq_id = self.trim_isoforms(seq_id,n=1)
            else:
                # There are no NPR-1-*
                # See if there is NPR-1. Look for that exact match
                matches = self.gpcr_names == name
                
                if np.sum(matches)>0:
                    seq_id = self.gpcr_seq_ids[matches][0]
                else:
                    # Nothing found. Try in the cache.
                    match = np.where(name.lower()==self.cached_seq_ids[:,0])[0]
                    if len(match)>0:
                        seq_id = self.cached_seq_ids[match[0],1]
                    else:
                        print("Need cached seq id for",name)
                        seq_id = name
            
            seq_ids.append(seq_id)
            
        if was_scalar: seq_ids = seq_ids[0]
        return seq_ids
   
    def get_gpcr_name_from_seq_id(self,seq_ids):
        if type(seq_ids) in [str,np.str_]: 
            seq_ids=[seq_ids]
            was_scalar = True
        else:
            was_scalar = False

        names = []
        for seq_id in seq_ids:
            # Check if seq_id is already trimmed
            seq_id_trimmed = len(seq_id.split("-"))==1
            seq_id_ = seq_id+"-"
            matches_ = np.flatnonzero(np.char.find(
                                            self.gpcr_seq_ids,seq_id_)!=-1)
            
            if len(matches_)>0:
                # There are NPR-1-1, etc. Use the matches of NPR-1-* BUT ALSO
                # the exact matches to NPR-1, which would otherwise be 
                # discarded.
                matches = matches_
                
                # No matter if there are NPR-1-1, if the name you passed has no
                # isoform label, then any of the NPR-1-* will give you the
                # sequence name you want, after trimming.
                name = self.gpcr_names[matches[0]]
            else:
                # There are no NPR-1-*
                # See if there is NPR-1. Look for that exact match
                matches = self.gpcr_seq_ids == seq_id
                
                if np.sum(matches)>0:
                    name = self.gpcr_seq_ids[matches][0]
                else:
                    # Nothing found. Try in the cache.
                    match = np.where(seq_id.lower()==self.cached_seq_ids[:,1])[0]
                    if len(match)>0:
                        name = self.cached_seq_ids[match[0],0]
                    else:
                        print("Need cached seq id for",seq_id)
                        name = ""
            
            if seq_id_trimmed:
                name = self.trim_isoforms(name)
            names.append(name)
            
        if was_scalar: names = names[0]
        return names
   
    @staticmethod     
    def trim_isoforms(b,n=2,delimiter="-",trim_letters=True):
        '''Trim isoform/variant label ("-n" or last letter) from a name.
        
        Parameters
        ----------
        b: str or array_like of str
            Name to be trimmed.
        n: int (optional)
            How many dash-separated parts of the string to expect in a 
            name without any isoform labeling (e.g. n=2 for FLP-1).
        delimiter: str (optional)
            Delimiter of the name sub-units. Default: "-"
        trim_letters: bool (optional)
            Whether to trim the last character of the name, if it is a letter.
            Default: True.
        
        Returns
        -------
        b_: str or array_like of str
            Trimmed name.
        '''
        if type(b) in [str,np.str_]:
            b = [b]
            was_scalar = True
        else:
            was_scalar = False
        
        b_ = []
        for a in b:
            l = len(a.split(delimiter))
            if l>n:
                b_.append(delimiter.join(a.split(delimiter)[:-(l-n)]))
            elif a[-1].isalpha():
                b_.append(a[:-1])
            else:
                b_.append(a)
        
        if was_scalar: b_ = b_[0]
        return b_
