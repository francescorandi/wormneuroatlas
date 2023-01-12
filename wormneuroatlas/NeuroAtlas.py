import numpy as np, matplotlib.pyplot as plt, json, h5py
import wormneuroatlas as wa

class NeuroAtlas:
    fname_neuron_ids = "neuron_ids.txt"
    '''Name of file containing full list of neurons'''
    fname_ganglia = "aconnectome_ids_ganglia.json"
    '''Name of the file containing the list of neurons divided in ganglia'''
    fname_signal_propagation = "funatlas.h5"
    '''Name of the file  containing the signal propagation atlas.'''
    fname_innexins = "GenesExpressing-unc-7-unc-9-inx-_-eat-5-thrs2.csv"
    '''Name of the file containing the expression levels of innexins'''
    fname_neuropeptides = "GenesExpressing-neuropeptides.csv"
    '''Name of the file containing the expression levels of neuropeptides'''
    fname_neuropeptide_receptors = "GenesExpressing-neuropeptide-receptors.csv"
    '''Name of the file containing the expression levels of neuropeptide ''' \
    '''receptors'''
    module_folder = ""
    '''Folder of the wormneuroatlas module'''
    
    aconn_sources = [
                     {"type": "whiteA", 
                     "fname": "aconnectome_white_1986_whole.csv"},
                     {"type": "whiteL4", 
                     "fname": "aconnectome_white_1986_L4.csv"},
                     {"type": "witvliet", 
                     "fname": "aconnectome_witvliet_2020_7.csv"},
                     {"type": "witvliet", 
                     "fname": "aconnectome_witvliet_2020_8.csv"}
                     ]
                     
    esconn_sources = [{"type": "bentley", 
                       "fname": "esconnectome_monoamines_Bentley_2016.csv", 
                       "transmitter_type": "monoamines"},
                      {"type": "bentley", 
                       "fname": "esconnectome_neuropeptides_Bentley_2016.csv", 
                       "transmitter_type": "neuropeptides"},
                     ]

    def __init__(self,merge_bilateral=False,merge_dorsoventral=False,
                 merge_numbered=False,merge_AWC=False,
                 load_connectomes=True, verbose=True,
                 *args,**kwargs):
        '''Class initialization.
        
        Parameters
        ----------
        merge_bilateral: bool (optional)
            Whether to merge bilateral pairs. Default: False.
        merge_dorsoventral: bool (optional)
            Whether to merge dorsoventral pairs. Default: False.
        merge_numbered: bool (optional)
            Whether to merge numbered neuron sets. Default: False.
        '''         
        
        if "\\" in wa.__file__:
            self.folder_sep = char = "\\"
        else:
            self.folder_sep = char = "/"
        self.module_folder = char.join(wa.__file__.split(char)[:-1])+char+\
                             "data"+char
        
        # Load the full list of neuronal ids from file       
        self._neuron_ids = np.loadtxt(self.module_folder+self.fname_neuron_ids,
                                     dtype=str)[:,1]
        
        # Set class options
        self.merge_bilateral = merge_bilateral  
        self.merge_dorsoventral = merge_dorsoventral
        self.merge_numbered = merge_numbered
        self.merge_AWC = merge_AWC
        
        if merge_numbered and verbose:
            print("Funatlas: Note that IL1 and IL2 will not be merged, as "+\
                   "well as VB1, VA1, VA11 because "+\
                   "they are different from the other VB and VA.")
        self.verbose = verbose
        
        # Compute the approximated (reduced) list of neuronal ids, e.g.
        # orphaning all the labels of their L/Rs.
        self.neuron_ids = self.approximate_ids(
                                self._neuron_ids,self.merge_bilateral,
                                self.merge_dorsoventral,self.merge_numbered,
                                self.merge_AWC)
        self.neuron_ids = np.unique(self.neuron_ids)
        self.n_neurons = self.neuron_n = len(self.neuron_ids)
        self.merger = self.get_merger()
        
        # Load classifications of neurons into head and pharynx, and 
        # repeat approximation done above.
        self.ganglia, self.head_ids, self.pharynx_ids = self.load_ganglia()
        self.head_ids = self.approximate_ids(
                                self.head_ids,self.merge_bilateral,
                                self.merge_dorsoventral,self.merge_numbered,
                                self.merge_AWC)
        self.head_ids = np.unique(self.head_ids)
        self.head_ai = self.ids_to_ai(self.head_ids)
        
        self.pharynx_ids = self.approximate_ids(
                                self.pharynx_ids,self.merge_bilateral,
                                self.merge_dorsoventral,self.merge_numbered,
                                self.merge_AWC)
        self.pharynx_ids = np.unique(self.pharynx_ids)
        self.pharynx_ai = self.ids_to_ai(self.pharynx_ids)
        
        if load_connectomes:
            self.load_aconnectome_from_file()
            self.load_extrasynaptic_connectome_from_file()
            self.aconn_loaded = True
            self.esconn_loaded = True
        else:
            self.aconn_loaded = False
            self.esconn_loaded = False
        self.load_signal_propagation_atlas()
        
        try:
            self.wormbase = wa.WormBase()
        except:
            # WormBase could not be loaded, likely because libcurl is missing.
            # Not printing anything here because the __init__.py does already.
            pass
        self.cengen = wa.Cengen()
        self.pepgpcr = wa.PeptideGPCR()
        
        ###########################################################
        # Manage conversions between different styles of neuron IDs
        ###########################################################
        
        # Using the function to convert CeNGEN-style neuron IDs to atlas-style
        # IDs, build the more useful reverse lookup to get the cengen_ai from 
        # the ai.
        # Get all the CeNGEN-style IDs
        cengen_ids = self.cengen.get_neuron_ids()
        # Convert them to all the corresponding atlas-style neuron IDs
        cengen_ids_to_atlas = self.cengen_ids_to_atlas(cengen_ids)
        # Initialize array of cengen_i given ai. In this way the CeNGEN 
        # matrices can be indexed with matrix[neura.cengen_i[ai]]
        self.cengen_is = -1*np.zeros(self.neuron_n,dtype=int)
        
        for ai in np.arange(self.n_neurons):
            for ci in np.arange(len(cengen_ids_to_atlas)):
                app_c_ids_to_atlas = self.approximate_ids(
                                cengen_ids_to_atlas[ci],self.merge_bilateral,
                                self.merge_dorsoventral,self.merge_numbered,
                                self.merge_AWC)
                if self.neuron_ids[ai] in app_c_ids_to_atlas:
                    self.cengen_is[ai] = ci
    
    
    def approximate_ids(self,ids,merge_bilateral=True,merge_dorsoventral=True,
                 merge_numbered=True,merge_AWC=False):
        '''Approximate the neural IDs, by orphaning them of the bilateral L/R
        and/or the dorso-ventral V/D and/or of the numbering for the 
        undistinguishable neurons in the retrovesicular ganglion and ventral
        nerve cord.
        
        Parameters
        ----------
        ids: str or array_like of str
            IDs to be approximated.
        merge_bilateral: bool (optional)
            Whether to merge left/right pairs. Default: True.
        merge_dorsoventral: bool (optional)
            Whether to merge dorsal/ventral pairs, triplets, and quadruples. 
            Default: True.
        merge_numbered: bool (optional)
            Whether to merge the undistinguishable merged neurons in the 
            retrovesicular ganglion and ventral nerve cord. Note that VB1, VA1, 
            VA11 will not be merged because they are different from the other 
            VB and VA. IL1 and IL2 will also not be merged. Default: True.
            
        Returns
        -------
        out_ids: str or list of str
            Approximated IDs. Single string if ids was a string, list of strings
            if ids was an array_like of strings.
        '''
        
        # If ids is a string, convert it to a list of string for processing.
        if type(ids)!=str:
            iter_input=True
        else: 
            ids = [ids] 
            iter_input=False
            
        out_ids = []
        for in_id in ids:
            if in_id is None:
                out_ids.append("")
                continue
            out_id = in_id
                    
            if len(in_id)>0:
                if in_id in ["AWCOFF","AWCON"]:
                    # Separate treatment: The anatomical data calls them AWCL 
                    # and AWCR, not ON and OFF. Special treatment also in
                    # self.ids_to_ai().
                    if merge_bilateral or merge_AWC:
                        out_id = "AWC"
                    #elif not merge_AWC and not merge_bilateral:
                    #    if in_id=="AWCOFF": out_id="AWCL"
                    #    elif in_id=="AWCON": out_id="AWCR"
                        
                if merge_bilateral and in_id not in ["AQR"]:
                    if in_id[-1]=="L":
                        if in_id[:-1]+"R" in self._neuron_ids or\
                         in_id[:-2]+"D"+"R" in self._neuron_ids or\
                         in_id[:-2]+"V"+"R" in self._neuron_ids or\
                         in_id[:-2]+"D" in self._neuron_ids: 
                            out_id = in_id[:-1]+"_"
                    elif in_id[-1]=="R":
                        if in_id[:-1]+"L" in self._neuron_ids or\
                         in_id[:-2]+"D"+"L" in self._neuron_ids or\
                         in_id[:-2]+"V"+"L" in self._neuron_ids or\
                         in_id[:-2]+"D" in self._neuron_ids:
                            out_id = in_id[:-1]+"_"
                            
                if merge_dorsoventral:
                    if len(out_id)==4:
                        if out_id[-1]=="D":
                            # To make SMBD -> SMB_ (because there is SMBVL)
                            
                            if out_id[:-1]+"V" in self._neuron_ids or\
                             out_id[:-1]+"VL" in self._neuron_ids or\
                             out_id[:-1]+"VR" in self._neuron_ids:
                                out_id = out_id[:-1]+"_"
                                
                        elif out_id[-1]=="V":
                            if out_id[:-1]+"D" in self._neuron_ids or\
                             out_id[:-1]+"DL" in self._neuron_ids or\
                             out_id[:-1]+"DR" in self._neuron_ids:
                                out_id = out_id[:-1]+"_"
                    
                    if len(out_id)==5:
                        if out_id[-2]=="V":
                            # To make OLQVL -> OLQ_L
                            if out_id[:-2]+"D"+out_id[-1] in self._neuron_ids or\
                              out_id[:-2]+"D"+"L" in self._neuron_ids or\
                              out_id[:-2]+"D"+"R" in self._neuron_ids or\
                              out_id[:-2]+"D" in self._neuron_ids:
                                out_id = out_id[:-2]+"_"+out_id[-1]
                                
                        elif out_id[-2]=="D":
                            # To make OLQVL -> OLQ_L
                            if out_id[:-2]+"V"+out_id[-1] in self._neuron_ids or\
                             out_id[:-2]+"V"+"L" in self._neuron_ids or\
                             out_id[:-2]+"V"+"R" in self._neuron_ids or\
                             out_id[:-2]+"V" in self._neuron_ids:
                                out_id = out_id[:-2]+"_"+out_id[-1]
                                
                            
                if merge_numbered:
                    if len(out_id)>2 and\
                      out_id not in ["VB1","VA1","VA11","IL1","IL2"]:
                        if out_id[-2:].isdigit():
                            out_id = out_id[:-2]+"_"
                        elif out_id[-1].isdigit():
                            out_id = out_id[:-1]+"_"
                            
                    
            out_ids.append(out_id)
            
        if iter_input: return out_ids
        else: return out_ids[0]
    
    def ids_to_ai(self,ids):
        '''Converts IDs to atlas indices.
        
        Parameters
        ----------
        ids: string or array_like of strings 
            IDs of the neurons.
        
        Returns
        -------
        atlas_i: int or numpy.ndarray of int
            Atlas-indices of the input IDs. Scalar or array depending on type of
            input.
        '''
        
        if type(ids)==str: ids = [ids]; was_string = True
        else: was_string = False
        
        n_ids = len(ids)
        
        atlas_i = -1*np.ones(n_ids,dtype=int)
        # Approximate the input identities according to the same settings as for
        # this atlas' instance.
        approx_ids = self.approximate_ids(
                            ids,self.merge_bilateral,self.merge_dorsoventral,
                            self.merge_numbered,self.merge_AWC)

        atlas_ids_B = np.copy(self.neuron_ids)
        for i in np.arange(len(atlas_ids_B)):
            atlas_ids_B[i] = atlas_ids_B[i].replace("_","")
            
        for i in np.arange(n_ids):
            approx_ids[i] = approx_ids[i].replace("_","")
                
        for i in np.arange(n_ids):
            matched = np.where(atlas_ids_B==approx_ids[i])[0]
            if matched.shape[0]>0:
                atlas_i[i] = matched[0]
        
        if was_string:
            return atlas_i[0]
        else:
            return atlas_i
            
    def ids_to_ais(self,*args,**kwargs):
        '''alias'''
        return self.ids_to_ai(*args,**kwargs)
        
    def get_merger(self):
        ''' Build the merger[i] = [k,l] that says that neurons k,l have
        been merged into neuron i.
        
        Returns
        -------
        merger: numpy.ndarray
            The length of each element of merger is the number of neurons that
            have been merged into each class.
        '''
        
        if self.merge_bilateral or self.merge_dorsoventral or \
           self.merge_numbered or self.merge_AWC:
            
            merger = np.empty(self.neuron_n,dtype=object)
            for i in np.arange(len(merger)): merger[i] = []
            
            for k in np.arange(len(self._neuron_ids)):
                ai = self.ids_to_ai([self._neuron_ids[k]])[0]
                merger[ai].append(k)
            
        else:
            merger = np.arange(self.neuron_n).reshape((self.neuron_n,1))
            
        return merger
        
    def merge(self, A, merger=None):
        '''Given an input array or matrix, the function merges it, based on the
        neuron-class approximation chosen at instantiation, by averaging the 
        elements of the matrices in the original reference frame. 
        
        Parameters
        ----------
        A: numpy.ndarray
            1-D or 2-D array to be merged. 
        
        merger: array_like
            merger[i] is the list [k,] of the neuron k that have been merged
            into class i. If None, A is merged based on the internally stored
            merger (that merges completely unmerged neurons into the approximate
            classes chosen at instantiation). Default: None.
            
        Returns
        -------
        B: array_like
            Merged version of A. Dimensionality and dtype are the same as for A.
        '''
        
        if merger is None:
            if not (self.merge_bilateral or self.merge_dorsoventral or\
                    self.merge_numbered or self.merge_AWC):
                # Nothing to be merged based on the approximations.
                return A
            else:
                merger = self.merger
        
        completeness = np.zeros(len(A),dtype=bool)
        if check_completeness:
            for m in merger:
                completeness[m] = True
        
        if not np.all(completeness):
            print("The merger is incomplete. Returning the unmerged array.")
            return A
        
        if len(A.shape)==1:
            B = np.zeros(len(merger),dtype=A.dtype)
            
            for i in np.arange(len(meger)):
                B[i] = np.average(A[np.array(merger[i])])
            
        elif len(A.shape)==2:
            B = np.zeros((len(merger),len(merger)),dtype=A.dtype)
            
            for i in np.arange(len(meger)):
                for j in np.arange(len(meger)):
                    b = B[np.array(merger[i])][:,np.array(merger[j])]
                    B[i,j] = np.average(b)
            
        else:
            raise ValueError("Input can only be 1 or 2 dimensional.")
            
        
        return B
            
            
    def load_ganglia(self):
        # Load the whole object
        g_f = open(self.module_folder+self.fname_ganglia,"r")
        ganglia = json.load(g_f)
        g_f.close()
        
        # Make a flattened list of neurons that are in the head
        head_ids = []
        for k in ganglia["head"]:
            for neu_id in ganglia[k]:
                head_ids.append(neu_id)
                
        pharynx_ids = []
        for k in ganglia["pharynx"]:
            for neu_id in ganglia[k]:
                pharynx_ids.append(neu_id)
        
        return ganglia, head_ids, pharynx_ids
    
    def reduce_to_head(self,A):
        '''Returns the submatrix or subarray of A corresponding to the head 
        neurons only.
        
        Parameters
        ----------
        A: numpy.ndarray
            Matrix or array to be cut.
        
        Returns
        -------
        A_new: numpy.ndarray
            Matrix reduced to the head indices.
        '''
        
        if len(A.shape)>=2:
            A_new = A[self.head_ai][:,self.head_ai]
        elif len(A.shape)==1:
            A_new = A[self.head_ai]
        return A_new 
            
    def reduce_to_pharynx(self,A):
        '''Returns the submatrix or subarray of A corresponding to the pharynx 
        neurons only.
        
        Parameters
        ----------
        A: numpy.ndarray
            Matrix or array to be cut.
        
        Returns
        -------
        A_new: numpy.ndarray
            Matrix reduced to the pharynx indices.
        '''
        
        if len(A.shape)==2:
            A_new = A[self.pharynx_ai][:,self.pharynx_ai]
        elif len(A.shape)==1:
            A_new = A[self.pharynx_ai]
        return A_new 
        
    def ai_to_head(self,ais):
        '''Translate whole-body atlas-indices to positions of those indices
        in the head-restricted reference frame.
        
        Parameters
        ----------
        ais: array_like of int
            Whole-body atlas-indices.
        
        Returns
        -------
        head_ais: numpy.ndarray of int
            Translated indices.        
        '''
        
        ais_head = np.zeros_like(ais)
        for i in np.arange(len(ais)):
            ai_head = np.where(self.head_ai==ais[i])[0]
            if len(ai_head)>0:
                ais_head[i] = ai_head[0]
            else:
                ais_head[i] = -1
                
        return ais_head
        
    def cengen_ids_to_atlas(self,ids):
        '''Convert neural IDs from the CeNGEN naming convention to the 
        naming convention used in this atlas.
        
        Parameters
        ----------
        ids: str or array_like of str
            Neural IDs in the CeNGEN convention.
        
        Returns
        -------
        names_out: str or list of str
            Neural IDs in this atlas' convention. string or list of string 
            depending on type of input.        
        '''
        
        if type(ids)==str:
            ids = [ids]
            was_scalar=True
        else:
            was_scalar=False
            
        
        names_out = []
        for iid in np.arange(len(ids)):
            if ids[iid]=="IL1":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["IL1_","IL1D_","IL1V_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["IL1"]
                else:
                    names = ["IL1R","IL1L","IL1DR","IL1DL","IL1VR","IL1VL"]
            elif ids[iid]=="IL2_DV":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["IL2D_","IL2V_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = []
                else:
                    names = ["IL2DR","IL2DL","IL2VR","IL2VL"]
            elif ids[iid]=="IL2_LR":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["IL2_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = []
                else:
                    names = ["IL2R","IL2L"]
            elif ids[iid]=="CEP":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["CEPD_","CEPV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["CEP__"]
                else:
                    names = ["CEPDR","CEPDL","CEPVR","CEPVL"]
            elif ids[iid]=="OLQ":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["OLQD_","OLQV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["OLQ__"]
                else:
                    names = ["OLQDR","OLQDL","OLQVR","OLQVL"]
            elif ids[iid]=="RMD_DV":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["RMDD_","RMDV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["RMD__"]
                else:
                    names = ["RMDDR","RMDDL","RMDVR","RMDVL"]
            elif ids[iid]=="RMD_LR":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["RMD_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["RMD_"]
                else:
                    names = ["RMDL","RMDR"]
            elif ids[iid]=="RME_DV":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["RMED","RMEV"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["RME_"]
                else:
                    names = ["RMED","RMEV"]
            elif ids[iid]=="RME_LR":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["RME_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["RME_"]
                else:
                    names = ["RMEL","RMER"]
            elif ids[iid]=="SMD":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SMDD_","SMDV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SMD__"]
                else:
                    names = ["SMDDR","SMDDL","SMDVR","SMDVL"]
            elif ids[iid]=="URY":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["URYD_","URYV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["URY__"]
                else:
                    names = ["URYDR","URYDL","URYVR","URYDL"]
            elif ids[iid]=="URA":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["URAD_","URAV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["URA__"]
                else:
                    names = ["URADR","URADL","URAVR","URAVL"]
            elif ids[iid]=="SAA":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SAAD_","SAAV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SAA__"]
                else:
                    names = ["SAADR","SAADL","SAAVR","SAAVL"]
            elif ids[iid]=="SAB":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SABD_","SABV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SAB__"]
                else:
                    names = ["SABDR","SABDL","SABVR","SABVL"]
            elif ids[iid]=="SIA":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SIAD_","SIAV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SIA__"]
                else:
                    names = ["SIADR","SIADL","SIAVR","SIAVL"]
            elif ids[iid]=="SIB":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SIBD_","SIBV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SIB__"]
                else:
                    names = ["SIBDR","SIBDL","SIBVR","SIBVL"]
            elif ids[iid]=="SMB":
                if self.merge_bilateral and not self.merge_dorsoventral:
                    names = ["SMBD_","SMBV_"]
                elif self.merge_bilateral and self.merge_dorsoventral:
                    names = ["SMB__"]
                else:
                    names = ["SMBDR","SMBDL","SMBVR","SMBVL"]
            elif ids[iid]=="VD_DD":
                if self.merge_numbered:
                    names = ["VD","DD"]
                else:
                    names = ["VD"+str(i) for i in np.arange(14)[1:]]
                    names += ["DD"+str(i) for i in np.arange(7)[1:]]
            elif ids[iid]=="AWC_ON":
                if self.merge_AWC:
                    names = ["AWC"]
                else:
                    names = ["AWCON"]
            elif ids[iid]=="AWC_OFF":
                if self.merge_AWC:
                    names = ["AWC"]
                else:
                    names = ["AWCOFF"]
            else:
                match=np.flatnonzero(np.char.find(self.neuron_ids,ids[iid])==0)
                names = self.neuron_ids[match].tolist()
                
            names_out.append(names)
        
        if was_scalar:
            return names_out[0]
        else:
            return names_out
            
    
    ##########################
    ##########################
    ##########################
    # SIGNAL PROPAGATION ATLAS
    ##########################
    ##########################
    ##########################
    
    def load_signal_propagation_atlas(self):
        self.funatlas_h5 = h5py.File(
                        self.module_folder+self.fname_signal_propagation,"r")
        
        neuron_ids = [n.decode("utf-8") for n in self.funatlas_h5["neuron_ids"][:]]
        
        for i in np.arange(len(self.neuron_ids)):
            if neuron_ids[i]!=self.neuron_ids[i] and neuron_ids[i]!="AWCOF":
                print("In NeuroAtlas signal propagation",
                      "Something wrong with neuron id",
                      neuron_ids[i],self.neuron_ids[i])
        
        # delta F/F averaged over a time window and trials
        self.dFF = {}
        self.dFF["wt"] = self.funatlas_h5["wt"]["dFF"][:]
        self.dFF["unc31"] = self.funatlas_h5["unc31"]["dFF"][:]
        self.dFF["wt"] = self.merge(self.dFF["wt"])
        self.dFF["unc31"] = self.merge(self.dFF["unc31"])
        
        # delta F/F averaged over a time window (each i,j contains all the 
        # corresponding trials)
        self.dFF_all = {}
        self.dFF_all["wt"] = self.funatlas_h5["wt"]["dFF_all"][:]
        self.dFF_all["unc31"] = self.funatlas_h5["unc31"]["dFF_all"][:]
        self.dFF_all["wt"] = self.merge(self.dFF_all["wt"])
        self.dFF_all["unc31"] = self.merge(self.dFF_all["unc31"])
        
        # q values 
        self.q = {}
        self.q["wt"] = self.funatlas_h5["wt"]["q"][:]
        self.q["unc31"] = self.funatlas_h5["unc31"]["q"][:]
        self.q["wt"] = self.merge(self.q["wt"])
        self.q["unc31"] = self.merge(self.q["unc31"])
        
        # equivalence-testing q values
        self.q_eq = {}
        self.q_eq["wt"] = self.funatlas_h5["wt"]["q_eq"][:]
        self.q_eq["unc31"] = self.funatlas_h5["unc31"]["q_eq"][:]
        self.q_eq["wt"] = self.merge(self.q_eq["wt"])
        self.q_eq["unc31"] = self.merge(self.q_eq["unc31"])
        if self.merge_bilateral or self.merge_dorsoventral or \
           self.merge_numbered or self.merge_AWC:
            print("IMPORTANT NOTE: Because of the neuron merging, q values",
                  "have been merged and AVERAGED. This is not mathematically",
                  "accurate from the statistical point of view. Do not use",
                  "q values like this. New implementation soon.")
        
        # occurrence matrix: each entry is the length of dFF_all[i,j]
        self.occ1 = {}
        self.occ1["wt"] = self.funatlas_h5["wt"]["occ1"][:]
        self.occ1["unc31"] = self.funatlas_h5["unc31"]["occ1"][:]
        self.occ1["wt"] = self.merge(self.occ1["wt"])
        self.occ1["unc31"] = self.merge(self.occ1["unc31"])
        
        # check that the ordering of the kernel keys is consistent with the 
        # definition used here
        if not self.funatlas_h5.attrs["kernels_keys"].decode("utf-8") == \
               "g,factor,power_t,branch":
            print("There is a problem with the ordering of the kernel keys.")
            print("DO NOT USE THE KERNELS BEFORE SOLVING THE PROBLEM.")
    
    def get_signal_propagation_map(self,strain="wt"):
        return self.dFF[strain]
        
    def get_sigpropmap(self,*args,**kwargs):
        '''alias
        '''
        return self.get_signal_propagation_map(*args,**kwargs)
    
    def get_signal_propagation_map_all(self,strain="wt"):
        return self.dFF_all[strain]
    
    def get_sigpropmap_all(self,*args,**kwargs):
        '''alias
        '''
        return self.get_signal_propagation_map_all(*args,**kwargs)
        
    def get_signal_propagation_occurrence_matrix(self,strain="wt"):
        return self.occ1[strain]
        
    def get_signal_propagation_q(self,strain="wt"):
        return self.q[strain]
        
    def get_kernel(self, i, j, strain="wt"):
        '''Return the average kernel for pair i<-j. Note that if neurons have 
        been merged upon instantiation of the NeuroAtlas class, the kernel
        that is returned is '''
        
        exp = []
        
        count = 0
        for ip in self.merger[i]:
            for jp in self.merger[j]:
                k = self.funatlas_h5[strain]["kernels"][i,j][:]
                if len(k)>0:
                    count += 1
                    k = k.reshape((len(k)//4,4))
                    
                    for k_ in k:
                        exp.append({"g":k_[0],"factor":k_[1],
                                     "power_t":k_[2],"branch":k_[3]})
                                     
        if len(exp>0):
            kernel = wa.ExponentialConvolution_min(exp)
            kernel.multiply_scalar_inplace(1.0/count)
            return kernel
        else:
            return None
        
            
    #######################
    #######################
    #######################
    # ANATOMICAL CONNECTOME
    #######################
    #######################
    #######################
        
    def load_aconnectome_from_file(self,chem_th=0,gap_th=0,exclude_white=False,
                                   average=False):
        self.aconn_chem, self.aconn_gap = \
                    self.get_aconnectome_from_file(chem_th,gap_th,exclude_white,
                                                   average=average)
        
    def get_aconnectome_from_file(self,chem_th=0,gap_th=0,exclude_white=False,
                                  average=False):
        '''Load the anatomical connectome data from all the sources listed in 
        the class.
        
        Parameters
        ----------
        chem_th: int (optional)
            Threshold on the counts of chemical synapses. Default: 0.
        gap_th: int (optional)
            Threshold on the counts of electrical synapses. Default: 0.
        exclude_white: bool (optional)
            Whether to exclude the older connectome from White et al. 1986.
            Default: False.
        average: bool (optional)
            Determines how the different connectomes are combined. If True,
            the average of the connectomes is taken, otherwise they are summed.
        
        
        Returns
        -------
        chem: numpy.ndarray
            chem[i,j] is the count of chemical synapses from j to i, averaged
            across the sources.
        gap: numpy.ndarray
            gap[i,j] is the count of gap junctions from j to i, averaged
            across the sources.
        '''
        chem = np.zeros((self.n_neurons, self.n_neurons))
        gap = np.zeros((self.n_neurons, self.n_neurons))
        
        # Reminder of the convention for AWCON and AWCOFF, since the anatomical
        # data has AWCL and AWCR.
        if not self.merge_AWC and not self.merge_bilateral and self.verbose:
            print("In loading the anatomical connectome, the following "+\
                  "conventions are used to allow for its comparison with the "+\
                  "other datasets: AWCL->AWCOFF and AWCR->AWCON")
        
        sources_used = 0
        for source in self.aconn_sources:
            if source["type"]=="white" and not exclude_white:
                c, g = self._get_aconnectome_white(
                                        self.module_folder+source["fname"],
                                        self.module_folder+source["ids_fname"])
            elif source["type"] in ["whiteL4","whiteA"] and not exclude_white:
                c, g = self._get_aconnectome_witvliet(
                                        self.module_folder+source["fname"])
            elif source["type"]=="witvliet":
                c, g = self._get_aconnectome_witvliet(
                                        self.module_folder+source["fname"])
            else:
                continue
            
            chem += c
            gap += g
            sources_used += 1
        
        if average:    
            chem /= sources_used
            gap /= sources_used
        
        chem[chem<=chem_th] = 0
        gap[gap<=gap_th] = 0
            
        return chem, gap
        
    def _get_aconnectome_white(self,fname,ids_fname):
        '''Load the anatomical connectome data from the file format for the
        White dataset.
        
        Returns
        -------
        chem: numpy.ndarray
            chem[i,j] is the count of chemical synapses from j to i.
        gap: numpy.ndarray
            gap[i,j] is the count of gap junctions from j to i.
        '''
        chem = np.zeros((self.n_neurons,self.n_neurons))
        gap = np.zeros((self.n_neurons,self.n_neurons))
        
        # Load the connectome matrices
        f = open(fname,"r")
        conn = json.load(f)
        f.close()
        
        # Load the associated ids
        white_ids = []
        f = open(ids_fname,'r')
        for l in f.readlines():
            neu_name = l.split("\t")[1]
            if neu_name[-1:] == "\n": neu_name = neu_name[:-1]
            white_ids.append(neu_name)
        f.close()
        
        a_is = self.ids_to_ai(white_ids)
        
        for i_w in np.arange(self.n_neurons):
            a_i = a_is[i_w]
            for j_w in np.arange(self.n_neurons):
                a_j = a_is[j_w]
                # The matrices have to be transposed, but they're lists
                # so switch i and j
                chem[a_i,a_j] += conn["chemical"][j_w][i_w]
                gap[a_i,a_j] += conn["electrical"][j_w][i_w]
                if conn["electrical"][i_w][j_w] == 0:
                    gap[a_j,a_i] += conn["electrical"][j_w][i_w]
                #chem[a_i,a_j] = conn["chemical"][i_w][j_w]
                #gap[a_i,a_j] = conn["electrical"][i_w][j_w]
                
        return chem, gap
        
    def _get_aconnectome_witvliet(self,fname):
        '''Load the anatomical connectome data from the file format for the
        Witvliet dataset, or datasets stored by Witvliet in Nemanode.
        
        Returns
        -------
        chem: numpy.ndarray
            chem[i,j] is the count of chemical synapses from j to i.
        gap: numpy.ndarray
            gap[i,j] is the count of gap junctions from j to i.
        '''
        chem = np.zeros((self.n_neurons,self.n_neurons))
        gap = np.zeros((self.n_neurons,self.n_neurons))
        
        f = open(fname,'r')
        lines = f.readlines()
        f.close()
        
        for l in lines[1:]:
            sl = l.split("\t")
            id_from = sl[0]
            id_to = sl[1]
            conn_type = sl[2]
            conn_n = int(sl[3])
            
            if conn_n!=0:
                ai_from = self.ids_to_ai(id_from)
                ai_to = self.ids_to_ai(id_to)
                
                if conn_type == "chemical":
                    chem[ai_to,ai_from] += conn_n
                elif conn_type == "electrical":
                    gap[ai_to,ai_from] += conn_n
                    #gap[ai_from,ai_to] += conn_n
                    
        for i in np.arange(self.n_neurons):
            for j in np.arange(self.n_neurons):
                gap[i,j] = max(gap[i,j],gap[j,i])
                
        return chem, gap
        
    def get_anatomical_connectome(self):
        return self.aconn_chem + self.aconn_gap
    
    def get_aconn(self,*args,**kwargs):
        return self.get_anatomical_connectome(*args,**kwargs)
        
    def get_gap_junctions(self):
        return self.aconn_gap
        
    def get_electrical_synapses(self,*args,**kwargs):
        return self.get_gap_junctions(*args,**kwargs)
        
    def get_chemical_synapses(self):
        return self.aconn_chem
        
    def get_effective_aconn(self,s=1):
        '''Returns the effective anatomical connectome computed with a given 
        gain for the connections.
        
        Parameters
        ----------
        s: float (optional)
            Gain, or modulation of the strength of the connections. Default: 1.
        
        Returns
        -------
        G: 2D numpy.ndarray
            Effective anatomical connectome.
        '''
        g = s*(self.aconn_chem+self.aconn_gap)
        diag_i = np.diag_indices(self.n_neurons)
        g[diag_i] = 0.0
        
        G = np.zeros_like(g)
        
        for j in np.arange(self.n_neurons):
            b = np.copy(g[:,j])
            a = np.copy(-g)
            a[:,j] = 0.0
            a[j,j] = -g[j,j]
            for i in np.arange(self.n_neurons): a[i,i] += 1.
            
            G[:,j] = np.linalg.solve(a,b)
            
        return G
        
    def get_effective_aconn2(self,max_hops=None):
        '''Returns the boolean anatomical connectome either at convergence or
        with a maximum number of hops.
        
        Parameters
        ----------
        max_hops: int (optional)
            Maximum number of hops to allow for the effective connectivity. If
            None, the effective connectivity is returned at convergence. The
            boolean single-hop anatomical connectome is returned with max_hops
            equal to 1. Default: None.
            
        Returns
        -------
        c: 2D numpy.ndarray of int
            Effective connectome. Values are either 0, if the neurons are in
            disjoint sets, or 1, if the neuron are connected.        
        '''
        c = (self.aconn_chem+self.aconn_gap)
        cl = np.min(c[c>0])
        c = np.clip(c,0,cl)/cl
        
        # The starting c is already with 1 hop (from 0 would be additional hops)
        p = 1
        if max_hops is None: max_hops = 10000
        while True:
            if p==max_hops: break
            old_c = c.copy()
            # c_ij = c_ik c_kj
            for i in np.arange(len(c)):
                for j in np.arange(len(c)):
                    c[i,j] = old_c[i,j]+np.sum(old_c[i,:]*old_c[:,j])
            c = np.clip(c,0,1)
            if np.all(c==old_c): break
            p +=1
            
        return c
        
    def get_effective_aconn3(self,max_hops=None,gain_1='average',c=None):
        '''Returns the anatomical connectome either at convergence or
        with a maximum number of hops, by setting the gain equal to 1 for 
        a given number of contacts.
        
        Parameters
        ----------
        max_hops: int (optional)
            Maximum number of hops to allow for the effective connectivity. If
            None, the effective connectivity is returned at convergence. The
            boolean single-hop anatomical connectome is returned with max_hops
            equal to 1. Default: None.
            
        gain_1: str or int (optional)
            Where to set the gain equal to 1 in terms of number of contacts.
            If it is set to 'average', the median number of contacts corresponds
            to gain=1. Default: average.
            
        Returns
        -------
        c: 2D numpy.ndarray of int
            Effective connectome.   
        '''
        if c is None:
            c = (self.aconn_chem+self.aconn_gap)
        
        if gain_1=='average':
            # Set the median number of connections to gain 1
            gain_1 = np.average(c)
            
        c /= gain_1
        orig_max = np.max(c)
            
        # The starting c is already with 1 hop (from 0 would be additional hops)
        p = 1
        if max_hops is None: max_hops = 10000
        while True:
            if p==max_hops: break
            old_c = c.copy()
            # c_ij = c_ik c_kj
            for i in np.arange(len(c)):
                for j in np.arange(len(c)):
                    c[i,j] = old_c[i,j]+\
                             np.sum(np.delete(old_c[i,:]*old_c[:,j],(i,j)))
            #c *= orig_max/np.max(c)
            if np.all(c==old_c): 
                print("Aconnectome converged in",p); break
            p +=1
            
        return c
        
    @staticmethod
    def _converged_matrix_power(A,n_max=100,eps=1e-2,return_all=False):
        A0 = A
        A_old = A
        n = 0
        e = 2*eps
        while True and n<n_max and e>eps:
            A = np.dot(A,A0)
            e = np.sum(np.absolute(A_old-A))/np.sum(A_old)
            A_old = A
            n += 1
        
        if return_all:
            return A,n
        else:
            return A
            
    @staticmethod
    def corr_from_eff_causal_conn(A):
        '''Built a correlation-type matrix from a causal-connectivity-type
        matrix.
        
        Parameters
        ----------
        A: 2D numpy.ndarary
            Causal-connectivity-type matrix.
            
        Returns
        -------
        B: 2D numpy.ndarray
            Correlation-type matrix.
        '''
        B = A+A.T
        C = np.zeros_like(B)
        count = np.zeros_like(B)
        
        n = A.shape[0]
        ns = np.arange(n)
        
        for i in np.arange(n):
            for k in np.arange(n):
                if k==i: continue
                for j in np.arange(n):
                    if j==i or j==k: continue
                    if np.isfinite(A[i,j]) and np.isfinite(A[k,j]):
                        C[i,k] += A[i,j]*A[k,j]
                        count[i,k] += 1
        B[count!=0] += C[count!=0]/count[count!=0]
        
        return B
    
    @staticmethod
    def symmetrize_nan_preserving(A):
        '''Symmetrize a matrix preserving the nan entries as such.
        '''
        nanmask = np.isnan(A)*np.isnan(A.T)
        A_ = 0.5*(np.nansum([A,A.T],axis=0))
        A_[nanmask] = np.nan
        
        return A_
            
    
    @staticmethod    
    def _get_next_hop_aconn(c):
        '''Compute the boolean map of strictly-2-hops connections given an input 
        connectome. To be used recursively in the calculation of the 
        strictly-n-hops connections.
        
        Parameters
        ----------
        c: numpy.ndarray of bool
            Input boolean connectome.
        
        Returns
        -------
        c_new: numpy.ndarray of bool
            Output strictly-2-hop boolean connectome.
        
        '''
        n = c.shape[0]
        
        c_new = np.zeros_like(c,dtype=bool)
        for j in np.arange(n):
            for i in np.arange(n):
                c_new[i,j] = not c[i,j] and np.any(c[i,:]*c[:,j])
                        
        return c_new
        
    def get_n_hops_aconn(self,n_hops):
        '''Return the anatomical connectome at n_hops from the original
        connectome.
        
        Parameters
        ----------
        n_hops: int
            Number of hops.
        
        Returns
        -------
        c_new: numpy.ndarray of bool
            Connectome n_hops away from the original one.
        '''
        
        c = self.get_boolean_aconn()
        if n_hops == 1: return c
        
        c_prev_hops = np.copy(c)
        for ih in np.arange(n_hops-1):
            # Get the connectome exactly 1-hop away from the current 
            # effective connectome.
            c_new = self._get_next_hop_aconn(c_prev_hops)
            # Update the current effective connectome by including all previous
            # connections and the new connections that you just found (any 
            # connection that has already been found corresponds to a connection
            # with less than the desired, exact n_hops).
            c_prev_hops = c_prev_hops+c_new
            
        return c_new
        
    def get_boolean_aconn(self):
        '''Returns the boolean anatomical connectome.
        
        Returns
        -------
        c: numpy.ndarray of bool
            Boolean anatomical connectome.
        '''
        
        c = (self.aconn_chem+self.aconn_gap) != 0
        return c
        
            
        
    def get_anatomical_paths(self,*args,**kwargs):
        '''Find all the paths between two neurons within a maximum numer of 
        hops.See Funatlas._get_paths() for the function arguments and returns. 
        In the  function arguments, skip the connectome matrix conn, which, 
        here, is set directly as the anatomical connectome.
        '''
        
        conn = self.aconn_chem + self.aconn_gap
        return self._get_paths(conn, *args, **kwargs)
                  
    def _get_paths(self,conn,i,j,max_n_hops=1,
                             return_ids=False,exclude_self_loops=True):
        '''Given a connectome, find all the paths between two neurons within a 
        maximum numer of hops.
        
        Parameters
        ----------
        conn: numpy.ndarray
            Connectome matrix.
        i: int or str
            Atlas-index or ID of the downstream neuron.
        j: int or str
            Atlas-index or ID of the upstream neuron.
        max_n_hops: int (optional)
            Maximum number of hops in which to perform the search. Default: 1.
        return_ids: bool (optional)
            Return also the paths with the IDs of the neurons, instead of their
            atlas-indices. Default: False.
        exclude_self_loops: bool (optional)
            Whether to exclude recurrent paths looping on the downstream
            neuron i.
            
        Returns
        -------
        paths: list of lists of integers
            paths[p][q] is the atlas index of q-th neuron on the p-th path. The 
            paths are ordered downstream-upstream.
        paths: list of lists of strings
            paths[p][q] is the index of q-th neuron on the p-th path. The 
            paths are ordered downstream-upstream. Returned only if return_ids
            is True.    
        '''
                             
        if type(i)==str: i = self.ids_to_ai(i)
        if type(j)==str: j = self.ids_to_ai(j)
        # You could add a check with a scalar Dyson equation to see if there
        # is a path at all. You can use the effective anatomical connectivity
        # at all steps to help speed up the process.
        paths_final = []
        
        paths_all = []
        # Populate paths_all with all the 1-hop connections that send signals
        # into i.
        for q in np.arange(self.n_neurons):
            if conn[i,q] != 0:
                paths_all.append([i,q])

        for h in np.arange(max_n_hops):
            paths_all_new = []
            for p in paths_all:
                if p[-1] == j: 
                    # Arrived at j
                    paths_final.append(p)
                elif h!=max_n_hops-1:
                    # Iterate over all the connections and add a hop
                    for q in np.arange(self.n_neurons):
                        if conn[p[-1],q]!=0 \
                            and not (exclude_self_loops and q==i):
                            new_p = p.copy()
                            new_p.append(q)
                            paths_all_new.append(new_p)
            paths_all = paths_all_new.copy()
        
        for p in paths_all:
            if p[-1] == j: 
                # Arrived at j
                paths_final.append(p)
                
        if return_ids:
            paths_final_ids = paths_final.copy()
            for i_p in np.arange(len(paths_final)):
                for q in np.arange(len(paths_final[i_p])):
                    paths_final_ids[i_p][q] = \
                        self.neuron_ids[paths_final_ids[i_p][q]]
            
            return paths_final, paths_final_ids
        else:
            return paths_final
            
    
    @staticmethod
    def _have_common_1_hop_upstream(c):
        '''Given a connectome, return the matrix that says whether two neurons
        have a common upstream neuron 1 hop away.
        
        Parameters
        ----------
        c: 2D numpy.ndarray
            Input connectome. Can be numeric or boolean.
        
        Returns
        -------
        c_new: 2D numpy.ndarray of bool
            c_new[i,j] is True if i and j have a common upstream neuron 1 hop 
            away in the original connectome c.
        '''
        n = len(c)
        c_new = np.zeros_like(c)
        
        for i in np.arange(n):
            for k in np.arange(n):
                # c_new[i,k] is True if there is at least one c[i,:] that is 
                # also in c[k,:].
                c_new[i,k] = np.sum(c[i,:]*c[k,:])!=0
                
        return c_new
        
    def have_common_n_hop_upstream(self,n=1):
        '''Return the matrix that says whether two neurons have a common 
        upstream neuron n hops away.
        
        Parameters
        ----------
        n: int (optional)
            Number of upstream hops at which to search for common upstream 
            neurons. Default 1.
        
        Returns
        -------
        c_new: 2D numpy.ndarray of bool
            c_new[i,j] is True if i and j have a common upstream neuron n hop 
            away in the original connectome c.
        '''
        if n!=1: raise ValueError("Implemented for n=1 only.")

        c = (self.aconn_chem+self.aconn_gap) != 0
        
        c_new = self._have_common_1_hop_upstream(c)
        
        return c_new
        
    def shuffle_aconnectome(self,shuffling_sorter=None):
        '''Shuffles the connectome internally stored in the instance of the 
        class.
        
        Parameters
        ----------
        shuffling_sorter: 1D numpy.ndarray (optional)
            Sorter that determines the suffling. Useful if the shuffler is 
            constructed somewhere else in the code. If None, a random shuffler 
            is generated. Default: None
        '''
        if shuffling_sorter is None:
            shuffling_sorter = self.get_shuffling_sorter()
            
        self.aconn_chem = self.shuffle_array(self.aconn_chem,shuffling_sorter)
        self.aconn_gap = self.shuffle_array(self.aconn_gap,shuffling_sorter)
        
    def load_innexin_expression_from_file(self):
        '''Load the expression of innexins from file (unc-7, unc-9, inx-*, 
        eat-5). Legacy, given the new interface to CeNGEN, but might be 
        still useful. 
        '''
        if not self.merge_bilateral:
            print("inx expression level available only with merge_bilateral")
            return None
        else:
            f = open(self.module_folder+self.fname_innexins,"r")
            lines = f.readlines()
            f.close()
            
            ids = lines[0][:-1].split(",")[4:]
            
            exp_levels = np.zeros((self.n_neurons,len(lines)-1))*np.nan
            genes = []
            
            for il in np.arange(len(lines)-1):
                l = lines[il+1]
                a=l.split(",")
                genes.append(a[1])
                
                for iid in np.arange(len(ids)):
                    names = self.cengen_ids_conversion([ids[iid]])[0]
                    for name in names:
                        aid = self.ids_to_ai(name)
                        if aid<0: continue
                        exp_levels[aid,il] = float(a[4+iid])
                    
            self.inx_exp_levels = exp_levels
            self.inx_genes = genes
                    
    def get_inx_exp_levels(self):
        unc7_i = self.inx_genes.index("unc-7")
        unc9_i = self.inx_genes.index("unc-9")

        unc7_9_to_others = np.zeros(self.n_neurons)
        
        for ai in np.arange(self.n_neurons):
            u79 = self.inx_exp_levels[ai,unc7_i]+self.inx_exp_levels[ai,unc9_i]
            if np.sum(self.inx_exp_levels[ai])>0:
                unc7_9_to_others[ai] = u79/np.sum(self.inx_exp_levels[ai])
        
        return self.inx_exp_levels, self.inx_genes, unc7_9_to_others
        
    def get_inx_exp_levels2(self,inx):
        '''Returns, for each neuron, the ratio of the expression level of a
        given innexin and the total expression of all the innexins.
        
        Parameters
        ----------
        inx: str
            Name of the innexin. Can be unc-7, unc-9, any inx-, or eat-5.
            
        Returns
        -------
        fr: 1D numpy.ndarray
            fr[ai] is the fractional expression of inx in neuron ai.
        '''
        
        inx_i = self.inx_genes.index(inx)
        
        fr = np.zeros(self.n_neurons)
        for ai in np.arange(self.n_neurons):
            tot_inx = np.sum(self.inx_exp_levels[ai])
            if tot_inx>0:
                fr[ai] = self.inx_exp_levels[ai,inx_i]/tot_inx
                
        return fr
        
    def get_fractional_gap_inx_mutants(self,mode='min'):
        '''Returns the remaining fraction of gap junctions, estimated by the
        product of the respective non-unc-7/9 fraction of innexins in the two
        neurons.        
        '''
        fr = np.ones((self.n_neurons,self.n_neurons))
        
        _,_,u79to = self.get_inx_exp_levels()
        
        if mode=="mult":
            fr *= (1-u79to[:,None])
            fr *= (1-u79to[None,:])
        else:
            for i in np.arange(fr.shape[0]):
                for j in np.arange(fr.shape[1]):
                    fr[i,j] = min( 1-u79to[i] , 1-u79to[j] )
        
        return fr
        
    def load_neuropeptide_expression_from_file(self):
        '''Legacy.'''
        if not self.merge_bilateral:
            print("neuropeptide expression level available only with",
                  "merge_bilateral.")
            return None
        else:
            f = open(self.module_folder+self.fname_neuropeptides,"r")
            lines = f.readlines()
            f.close()
            
            ids = lines[0][:-1].split(",")[4:]
            
            exp_levels = np.zeros((self.n_neurons,len(lines)-1))*np.nan
            genes = []
            
            for il in np.arange(len(lines)-1):
                l = lines[il+1]
                a=l.split(",")
                genes.append(a[1])
                
                for iid in np.arange(len(ids)):
                    names = self.cengen_ids_conversion([ids[iid]])[0]
                    
                    for name in names:
                        aid = self.ids_to_ai(name)
                        if aid<0: continue
                        if np.isnan(exp_levels[aid,il]):
                            exp_levels[aid,il] = float(a[4+iid])
                        else:
                            exp_levels[aid,il] += float(a[4+iid])
                    
            self.npt_exp_levels = exp_levels
            self.npt_genes = genes
            
    def load_neuropeptide_receptor_expression_from_file(self):
        '''Legacy.'''
        if not self.merge_bilateral:
            print("neuropeptide receptor expression level available only with",
                  "merge_bilateral")
            return None
        else:
            f = open(self.module_folder+self.fname_neuropeptide_receptors,"r")
            lines = f.readlines()
            f.close()
            
            ids = lines[0][:-1].split(",")[4:]
            
            exp_levels = np.zeros((self.n_neurons,len(lines)-1))*np.nan
            genes = []
            
            for il in np.arange(len(lines)-1):
                l = lines[il+1]
                a=l.split(",")
                genes.append(a[1])
                
                for iid in np.arange(len(ids)):
                    names = self.cengen_ids_conversion([ids[iid]])[0]
                    
                    for name in names:
                        aid = self.ids_to_ai(name)
                        if aid<0: continue
                        if np.isnan(exp_levels[aid,il]):
                            exp_levels[aid,il] = float(a[4+iid])
                        else:
                            exp_levels[aid,il] += float(a[4+iid])
                    
            self.nptr_exp_levels = exp_levels
            self.nptr_genes = genes
    
    ###########################################################
    ###########################################################
    ###########################################################
    # INTERFACE TO CENGEN AND NEUROPEPTIDE/GPCR DEORPHANIZATION
    ###########################################################
    ###########################################################
    ###########################################################
    
    def get_gene_expression(self, neuron_ais=None, neuron_ids=None, 
                            *args, **kwargs):
        '''Get gene expression. Relays the function call to 
        Cengen.get_gene_expression() after converting the neuron atlas-indices
        or neuron IDs to the corresponding CeNGEN-style index or ID.
        
        Parameters
        ----------
        *args, **kwargs:
            See Cengen.get_gene_expression()
        
        Returns
        -------
        gene_exp: 2D numpy.ndarray
            gene_exp[i,k] is the gene expression of the k-th gene in the i-th
            neuron.
        
        '''
        # Generate the list of CeNGEN indices (cis) of the requeste neurons.
        if neuron_ais is not None:
            neuron_cis = self.cengen_is[neuron_ais]
        elif neuron_ids is not None:
            try: neuron_ids[0] 
            except: neuron_ids = np.array([neuron_ids])
            neuron_ais = self.ids_to_ais(neuron_ids)
            # Check that all neuron IDs were found with the current IDs 
            # approximation.
            for q in np.arange(len(neuron_ais)):
                if neuron_ais[q]==-1:
                    e=neuron_ids[q]+" is not a valid ID with the current"+\
                    "neuron IDs approximation:\n"+\
                    "merge_bilateral "+str(self.merge_bilateral)+"\n"+\
                    "merge_dorsoventral "+str(self.merge_dorsoventral)+"\n"+\
                    "merge_numbered "+str(self.merge_numbered)+"\n"+\
                    "merge_AWC "+str(self.merge_AWC)+"\n"+\
                    "If you did not set any neuron IDs approximation, "+\
                    "you need to use individual neuron IDs, like AVAL "+\
                    "and not AVA."
                    raise ValueError(e)
                    
            neuron_cis = self.cengen_is[neuron_ais]
        else:
            neuron_cis = self.cengen_is
            
        kwargs["neuron_cis"] = neuron_cis
        
        gene_exp = self.cengen.get_expression(*args,**kwargs)
        
        return gene_exp
        
    def supplement_peptide_gpcr(self):
        self.pepgpcr.supplement()
        print("Supplementing peptide/GPCR combinations with older datasets.")
        
        return None
        
    def get_gpcrs_binding_to(self,*args,**kwargs):
        return self.pepgpcr.get_gpcrs_binding_to(*args,**kwargs)
        
    def get_peptides_binding_to(self,*args,**kwargs):
        return self.pepgpcr.get_peptides_binding_to(*args,**kwargs)
        
    def get_gpcr_names(self,*args,**kwargs):
        return self.pepgpcr.get_gpcr_names(*args,**kwargs)
        
    def get_peptide_names(self,*args,**kwargs):
        return self.pepgpcr.get_peptide_names(*args,**kwargs)
    
    def get_gpcr_seq_id_from_name(self,*args,**kwargs):
        return self.pepgpcr.get_gpcr_seq_id_from_name(*args,**kwargs)
    
    def get_gpcr_name_from_seq_id(self,*args,**kwargs):
        return self.pepgpcr.get_gpcr_name_from_seq_id(*args,**kwargs)
        
    def get_peptidergic_connectome(self,neuron_ids_from=None,neuron_ids_to=None,
                                   return_combos=False):
        '''Returns the neuronpeptidergic connectome based on CeNGEN peptide and
        GPCR expression data and the peptide/GPCR binding screen. 
        
        Parameters
        ----------
        neuron_ids_from: array_like of str or None (optional)
            Array of the IDs of the upstream neuron from which neuropeptide are
            released. If None, all neurons are included. Default: None.
        neuron_ids_to: array_like of str or None (optional)
            Array of the IDs of the downstream neuron to which the neuropeptide
            signal. If None, all neurons are includede. Default: None.
        return_combo: bool (optional)
            Whether to return also the list of peptide/GPCR combinations for
            each neuron pair. Default: False.
            
        Returns
        -------
        pep_esconn: 2D numpy.ndarray of int
            Peptidergic extrasynaptic connectome. pep_esconn[i,j] is the number
            of peptide/GPCR combinations for the connection i<-j. Note that i
            and j are indices in the neuron_ids_from and neuron_ids_to arrays,
            so if either of those arguments is not None, i or j might not 
            correspond to NeuroAtlas indices.
        pep_gpcr_combo: 2D numpy.ndarray of lists
            Peptide/GPCR combinations. pep_gpcr_combo[i,j] is a list of list
            containing the peptide/GPCR combinations for the connection i<-j.
            For example, pep_gpcr_combo[i,j] = [['FLP-1','FLPR-1'],[...,...]].
            Only returned if return_combo is True.        
        '''
        
        if neuron_ids_to is None:
            neuron_ids_to = self.neuron_ids
        if neuron_ids_from is None:
            neuron_ids_from = self.neuron_ids
        
        n_to = len(neuron_ids_to)
        n_from = len(neuron_ids_from)
        
        # Get the list of all peptide and GPCR names contained in the 
        # deorphanization screen.
        peptide_names = np.array(self.get_peptide_names(trim_isoforms=True))
        gpcr_names = np.array(self.get_gpcr_names(trim_isoforms=True))
        gpcr_seq_ids = np.array(self.get_gpcr_seq_id_from_name(gpcr_names))
                                                    
        # Get the neuron-resolved expression of those peptides and GPCR. For the
        # GPCRs, supplement the names with the sequence IDs, as some entries use
        # the sequence ID instead of the name.
        pep_exp = self.get_gene_expression(gene_names=peptide_names,
                                           neuron_ids=neuron_ids_from)
        gpcr_exp = self.get_gene_expression(gene_names=gpcr_names,
                                            gene_seq_ids=gpcr_seq_ids,
                                            neuron_ids=neuron_ids_to)
                                            
        # Make the extrasynaptic connectome (esconn) map for neuropeptide/GPCR
        pep_esconn = np.zeros((n_to,n_from))
        # Store also the peptide/GPCR combinations.
        if return_combos:
            pep_gpcr_combo = np.empty((n_to,n_from),dtype=object)
        
        for j in np.arange(n_from):
            # Find the peptides expressed in the upstream neuron j
            pep_expressed = peptide_names[pep_exp[j]!=0]
            # Find the GPCRs that bind to those peptides. Get the sequence IDs 
            # instead of the names.
            gpcr_binding_seq_ids = \
                  self.get_gpcrs_binding_to(pep_expressed,return_seq_ids=True)
            
            for i in np.arange(n_to):
                # Find the GPCRs expressed in the downstream neuron i
                gpcr_expressed = gpcr_seq_ids[gpcr_exp[i]!=0]
                
                # Count the number of GPCRs expressed in the downstream neuron 
                # that are also receptors for the peptides expressed in the 
                # upstream neuron. gpcrs requires two nested iterations because 
                # gpcrs[neuropeptide] is a list of all gpcrs binding to that 
                # neuropeptide.
                if return_combos: pep_gpcr_combo[i,j] = []
                for pep_i in np.arange(len(gpcr_binding_seq_ids)):
                    for gp_i in np.arange(len(gpcr_binding_seq_ids[pep_i])):
                        gp = gpcr_binding_seq_ids[pep_i][gp_i]
                        if gp in gpcr_expressed: 
                            pep_esconn[i,j]+=1
                            if return_combos:
                                pep_gpcr_combo[i,j].append(
                                                [pep_expressed[pep_i],gp])
                            
        if return_combos:
            return pep_esconn, pep_gpcr_combo
        else:
            return pep_esconn
        
        
    
    ##############################################
    ##############################################
    ##############################################
    # EXTRASYNAPTIC CONNECTOME FROM BENTLEY ET AL.
    ##############################################
    ##############################################
    ##############################################
    
    def load_extrasynaptic_connectome_from_file(self, *args, **kwargs):
         esconn_ma, esconn_np = self.get_extrasynaptic_connectome_from_file(
                                                *args,**kwargs)
         self.esconn_ma = esconn_ma
         self.esconn_np = esconn_np
    
    def get_extrasynaptic_connectome_from_file(self, transmitter_types=None,
                                               *args,**kwargs):
        '''Load the boolean extrasynaptic connectome data from all the sources 
        listed in the class.
        
        Parameters
        ----------
        transmitter_types: str or list of str (optional)
            Types of transmitters to allow (monoamines, neuropeptides, ...).
            If None, all are selected. Default: None.
        args, kwargs
        
        Returns
        -------
        ma: numpy.ndarray of bool
            Monoamine extrasynaptic connectome.
        np: numpy.ndarray of bool
            Neuropeptide extrasynaptic connectome.
        
        '''
        esconn_ma = np.zeros((self.n_neurons, self.n_neurons),dtype=bool)
        esconn_np = np.zeros((self.n_neurons, self.n_neurons),dtype=bool)
        
        if transmitter_types is not None:
            if type(transmitter_types)==str:
                transmitter_types = [transmitter_types]
        
        sources_used = 0
        for source in self.esconn_sources:
            if transmitter_types is not None:
                if source["transmitter_type"] not in transmitter_types:
                    continue
            if source["type"]=="bentley":
                esc = self._get_esconnectome_bentley(
                            self.module_folder+source["fname"],
                            *args,**kwargs)
            else:
                continue
            
            if source["transmitter_type"] == "monoamines":
                esconn_ma = np.logical_or(esconn_ma,esc)
            elif source["transmitter_type"] == "neuropeptides":
                esconn_np = np.logical_or(esconn_np,esc)
            
            sources_used += 1
        
        return esconn_ma, esconn_np
    
    def _get_esconnectome_bentley(self,fname,transmitters=None,receptors=None):
        '''Returns the extrasynaptic connectome from Bentley et al. 2016 "The 
        multilayer connectome of Caenorhabditis elegans" PLOS Comp. Bio.
        
        Parameters
        ----------
        fname: str
            Name of the csv file.
        transmitter: str or list of str (optional)
            Requested trasmitters. If None, no restriction is applied. 
            Default: None.
        receptors: str or list of str (optional)
            Requested receptors. If None, no restriction is applied. 
            Default: None.
        
        Returns
        -------
        esc: numpy.ndarray of bool
            Extrasynaptic connectome given the requested transmitters and
            receptors.
            
        '''
        esc = np.zeros((self.n_neurons,self.n_neurons),dtype=bool)
        
        if transmitters is not None:
            if type(transmitters)==str:transmitters = [transmitters]
        if receptors is not None:
            if type(receptors)==str:receptors = [receptors]
                
        
        f = open(fname,'r')
        lines = f.readlines()
        f.close()
        
        for l in lines[1:]:
            sl = l.split(",")
            id_from = sl[0]
            id_to = sl[1]
            trans = sl[2]
            recept = sl[3]
            if recept[-1] == "\n": recept = recept[-1]
            
            # Skip the line if the transmitter/receptor are not the requested
            # ones.
            if transmitters is not None:
                if trans not in transmitters: continue
            if receptors is not None:
                if recept not in receptors: continue
            
            ai_from = self.ids_to_ai(id_from)
            ai_to = self.ids_to_ai(id_to)
            
            esc[ai_to,ai_from] = True
                
        return esc
    
    def have_common_n_hop_es_upstream(self,n=1):
        if n!=1: raise ValueError("Implemented for n=1 only.")

        c = self.get_esconn()
        
        c_new = self._have_common_1_hop_upstream(c)
        
        return c_new
    
    
    def get_esconn(self):
        print("get_esconn changed ^ to +")
        return self.esconn_ma+self.esconn_np
        
    def get_monoaminergic_connectome(self):
        return self.esconn_ma
        
    def get_effective_esconn(self,maxit=100):
        escon = self.get_esconn()
        eescon = np.copy(escon)
        old_eescon = np.copy(escon)
        it = 0
        while True and it<maxit:
            for j in np.arange(escon.shape[0]):
                for i in np.arange(escon.shape[0]):
                    if escon[i,j]:
                        eescon[escon[:,i],j] = True
            if np.all(old_eescon==eescon): break
            old_eescon = eescon
            it+=1
        return eescon
        
    def get_extrasynaptic_paths(self, *args, **kwargs):
        '''Returns the extrasynaptic paths between neurons. See 
        Funatlas._get_paths() for the function arguments and returns. In the 
        function arguments, skip the connectome matrix conn, which, here,
        is set directly as the extrasynaptic connectome.
        '''
        
        conn = self.get_esconn()
        return self._get_paths(conn, *args, **kwargs)
        
    def shuffle_esconnectome(self,shuffling_sorter=None):
        if shuffling_sorter is None:
            shuffling_sorter = self.get_shuffling_sorter()
            
        self.esconn_ma = self.shuffle_array(self.esconn_ma,shuffling_sorter)
        self.esconn_np = self.shuffle_array(self.esconn_np,shuffling_sorter)
        
    ##########################
    ##########################
    ##########################
    # UTILITIES
    ##########################
    ##########################
    ##########################   
    
    @staticmethod
    def threshold_to_sparseness(A,sparseness,absolute=True,max_i=100):
        '''Returns the number by which a given matrix or tensor must be 
        thresholded to make it sparse as specified. Note, the tensor is not
        sparsified in place by this function.
        
        Parameters
        ----------
        A: array_like
            Tensor to be sparsified.
        sparseness: float
            Desired sparseness, expressed as fraction of non-zero elements.
        absolute: bool (optioal)
            Whether to sparsify and threshold based on absolute values.
            Default: True.
        max_i: int (optional)
            Maximum number of iterations allowed in estimating the thresold.
            Default: 100.
        
        Returns
        -------
        th: float
            Threshold that will sparsify the tensor A to the specified 
            sparseness.
        '''
        amax = np.max(A)
        amin = np.min(A)
        
        a0 = amin
        a1 = amax
        
        if absolute:
            B = np.absolute(A)
        else:
            B = A
            
        totB = np.prod(B.shape)
        
        i = 0
        while True:
            th = 0.5*(a0+a1)
            sp = np.sum(B>th)/totB
            if np.absolute(sp-sparseness)/sparseness < 1e-2 or i>max_i:
                break
            
            if sp<sparseness:
                a1 = th
            else:
                a0 = th
            i+=1
            
        return th
        
    ##########################
    ##########################
    ##########################
    # PLOTTING
    ##########################
    ##########################
    ##########################   
    
    @staticmethod
    def make_alphacolorbar(cax,vmin,vmax,tickstep,
                       alphamin,alphamax,nalphaticks,
                       cmap="viridis",bg_gray=0.2,around=0,lbl_lg=False,
                       lbl_g=False,alphaticklabels=None):
        '''Make a 2D colorbar with transparency on one axis.
        
        Parameters
        ----------
        cax: matplotlib axis
            Axis in which to draw the colorbar.
        vmin, vmax: float
            Minimum and maximum values of the colormap axis (y).
        tickstep: float
            Step of the ticks on the colormap axis (y).
        alphamin, alphamax: float
            Minimum and maximum alphas.
        nalphaticks: int
            Number of ticks on the alpha axis (x).
        cmap: matplotlib-accepted colormap (optional)
            Colormap.
        around: int (optional)
            np.around -ing precision for ticks. Default: 0.
        '''
        
        nticks = int((vmax-vmin)/tickstep)
        
        alphacm1 = np.zeros((200,200))
        alphacm1[:] = np.arange(alphacm1.shape[0])[:,None]
        alphacm2 = np.zeros_like(alphacm1)
        alphacm2[:] = (np.arange(alphacm1.shape[1])/alphacm1.shape[1])[None,:]
        
        background = np.ones_like(alphacm1)*bg_gray
        
        cax.yaxis.tick_right()
        cax.yaxis.set_label_position("right")
        cax.imshow(background,cmap="Greys",vmin=0,vmax=1)
        cax.imshow(alphacm1,alpha=alphacm2,aspect="auto",cmap=cmap,
                   interpolation="nearest",origin="lower")

        cax.set_xticks([0,alphacm1.shape[1]])
        if alphaticklabels is None:
            cax.set_xticklabels([str(alphamin),str(alphamax)])
        else:
            cax.set_xticklabels([alphaticklabels[0],alphaticklabels[1]])
        yticks = np.linspace(0,alphacm1.shape[0],nticks+1)
        cax.set_yticks(yticks)
        d = (vmax-vmin)/nticks
        yticklbl = [str(np.around(i*d+vmin,around)) \
                       for i in np.arange(len(yticks))]
        if lbl_lg: 
            yticklbl[0] = "<"+yticklbl[0]
            yticklbl[-1] = ">"+yticklbl[-1]
        if lbl_g:
            yticklbl[-1] = ">"+yticklbl[-1]
        cax.set_yticklabels(yticklbl)
        cax.set_xlim(0,alphacm1.shape[1])
        cax.set_ylim(0,alphacm1.shape[0])
        
    @classmethod
    def plot_matrix(cls, A, ids, alphas=None, 
                    cmap="Spectral_r", vmax=None, vmin=None,
                    black_diag=True, 
                    fig=None, figsize=(20,20), 
                    label="", alphalabel="",
                    labelx="stimulated", labely="responding",
                    labelsize=30,ticklabelsize=6,
                    alphamin=0.0,alphamax=1.0):
        '''Plot a matrix, add IDs as tick labels, and optionally use a second
        matrix to set gray levels to indicate, for example, confidence on the
        values in the first matrix. If this second matrix is passed, the 
        function generates a 2D colorbar that fades into gray for low 
        confidence.
        
        Parameters
        ----------
        A: 2D numpy.ndarray
            Matrix to be plotted.
        ids: array_like of str
            IDs to use as tick labels. They are used for both the x and y axes.
        alphas: 2D numpy.ndarray (optional)
            If not None, it is used to set the gray level of each entry in the
            matrix A. For example, it can represent the statistical confidence
            of the values in the matrix A. Default: None.
        cmap: str (optional)
            Colormap to use. Default: Spectral_r
        black_diag: bool (optional)
            Whether to black out diagonal entries. Default: True.
        fig: matplotlib figure (optional)
            If not None, this figure will be used for plotting. Default: None.
        figsize: array_like (2,) (optional)
            Figure size. Default: (20,20)
        label: str (optional)
            Label that will appear next to the colorbar. Default: "".
        alphalabel: str (optional)
            Label of the alpha/gray level. Default: "".
        labelx: str (optional)
            Label of the x axis. Default: "Stimulated".
        labely: str (optional)
            Label of the y axis. Default: "Responding".
        labelsize: int (optional)
            Label size. Default: 30.
        ticklabelsize: int (optional)
            Size of the tick labels. Default: 6.
        alphamin: float (optional)
            Minimum value for the clipping of the alpha values. Default: 0.0.
        alphamax: float (optional)
            Maximum value for the clipping of the alpha values. Default: 1.0.
            
        Returns
        -------
        fig: matplotlib figure
            Matplotlib figure containing the plot.
        ax: matplotlib axis
            Axis containing the plot.
        cax: matplotlib axis
            Axis containing the 2D colorbar. Returned only if alphas is not 
            None.
        '''            
        
                    
        if fig is None:
            cfn = plt.gcf().number
            if len(plt.gcf().axes) != 0: cfn += 1
            fig = plt.figure(cfn,figsize=figsize)
            
        if alphas is None:
            # There is no need for gridspec
            ax = fig.add_subplot(111)
        else:
            # Allow for more space 
            gs = fig.add_gridspec(1,10)
            ax = fig.add_subplot(gs[0,:9])
            cax = fig.add_subplot(gs[0,9:])
            
            # Adjust based on alphamin and alphamax
            alphas = np.clip((alphas-alphamin)/(alphamax-alphamin),0,1)
            
        if vmax is None:
            vmax = np.nanmax(A)
        
        if vmin is None:
            vmin = -vmax
        
        ax.imshow(0.*np.ones_like(A),cmap="Greys",vmax=1,vmin=0)
        blank_A = np.copy(A)
        blank_A[~np.isnan(A)] = 0.1
        ax.imshow(blank_A,cmap="Greys",vmin=0,vmax=1)
            
        im = ax.imshow(A,cmap=cmap,vmin=vmin,vmax=vmax,
                       alpha=alphas,interpolation="nearest")
        
        if black_diag:
            diagonal = np.diag(np.diag(np.ones_like(A)))
            new_diagonal = np.zeros_like(A)
            new_diagonal[np.where(diagonal == 1)] = 1
            ax.imshow(new_diagonal, cmap="binary", vmin=0, vmax=1, 
                      alpha=new_diagonal, interpolation="nearest")
            
        if alphas is None:
            plt.colorbar(im, label=label)
        else:
            cls.make_alphacolorbar(cax,vmin,vmax,0.1,alphamin,alphamax,2,
                                    cmap=cmap,around=1,lbl_lg=True)
            cax.set_xlabel(alphalabel,fontsize=labelsize//2)
            cax.set_ylabel(label,fontsize=labelsize//2)
            cax.tick_params(labelsize=labelsize//2)
        
        ax.set_xlabel(labelx,fontsize=labelsize)
        ax.set_ylabel(labely,fontsize=labelsize)
        
        ax.set_xticks(np.arange(len(ids)))
        ax.set_yticks(np.arange(len(ids)))
        ax.set_xticklabels(ids,fontsize=ticklabelsize,rotation=90)
        ax.set_yticklabels(ids,fontsize=ticklabelsize)
        
        ax.tick_params(axis="x", bottom=True, top=True, 
                       labelbottom=True, labeltop=True)
        ax.tick_params(axis="y", left=True, right=True, 
                       labelleft=True, labelright=True)
        ax.set_xlim(-0.5,len(ids)+0.5)
        
        if alphas is not None:
            return fig,ax,cax
        else:
            return fig,ax
        

