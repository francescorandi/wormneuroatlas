import json, numpy as np, urllib3
import wormneuroatlas as wa

class WormBase:
    
    gene_wbids_fname = "c_elegans.PRJNA13758.WS286.geneIDs.txt"
    module_folder = ""
    '''Folder of the wormneuroatlas module'''
    
    def __init__(self, gene_subset = None):
        if "\\" in wa.__file__:
            self.folder_sep = char = "\\"
        else:
            self.folder_sep = char = "/"
        self.module_folder = char.join(wa.__file__.split(char)[:-1])+char+\
                             "data"+char
        
        self.load_gene_wbids()
        
        if gene_subset is not None:
            # Make a gene_wbid_subset
            self.gene_wbid_subset = np.zeros(len(gene_subset))
            for i in np.arange(len(gene_subset)):
                if type(gene_subset[i]) != int:
                    if gene_subset[0][0:6] == "WBGene":
                        gene = self.wbid_to_wbint(gene_subset[i])
                    else:
                        gene = self.name_to_gene_wbid(gene_subset[i],alt=True)
                
    def load_gene_wbids(self):
        '''Load gene IDs, names, other codes, and types from file. See 
        WormBase.gene_wbids_fname for which version of the database is being 
        used.
        '''
        f = open(self.module_folder+self.gene_wbids_fname,"r")
        lines = f.readlines()
        f.close()
        self.g_wbids = []
        self.g_names = []
        self.g_seq_id = []
        self.g_types = []
        for l in lines:
            la = l.split(",")
            if la[4] != "Live": continue
            self.g_wbids.append(la[1])
            self.g_names.append(la[2])
            self.g_seq_id.append(la[3])
            self.g_types.append(la[5])
            
    def alt_name_to_gene_wbid(self,alt_name):
        '''Convert alternative gene names to gene ID.
        
        Parameters
        ----------
        alt_name: str
            Alternative name of gene.
            
        Returns
        -------
        g_wbid: int
            WormBase gene ID corresponding to the alternative gene name.
        '''
        
        try: self.alt_name_gene_wbid_warning
        except: 
            self.alt_name_gene_wbid_warning = True
            print("WormAtlas alt_name_to_gene_wbid() to be implemented. None.")
        g_wbid = None
        
        return g_wbid
    
    def name_to_gene_wbid(self,name,alt=False,dtype=int):
        '''Returns WormBase gene ID given either the official name or 
        alternative names. Specify alternative=True if the name is an 
        alternative name, or if you don't know. The function might run more 
        slowly if alternative==True.
        
        Parameters
        ----------
        name: str
            Gene name.
        alt: bool (optional)
            Whether to also consider alternative namings.
        
        Returns
        -------
        g_wbid: int
            WormBase gene ID corresponding to the gene name. g_wbid is None
            if no match has been found.
        '''
        try: 
            g_wbid_i = self.g_names.index(name)
            g_wbid = self.g_wbids[g_wbid_i]
        except: g_wbid = None
        
        if g_wbid is None and alt:
            g_wbid = self.alt_name_to_gene_wbid(name)
            
        if dtype==int and g_wbid is not None:
            g_wbid = self.wbid_to_wbint(g_wbid)
            
        return g_wbid
    
    @staticmethod
    def wbid_to_wbstr(gene_wbid):
        return "WBGene"+"{:>08}".format(str(gene_wbid))
        
    @staticmethod
    def wbid_to_wbint(gene_wbid):
        return int(gene_wbid[6:])
    
    @staticmethod
    def http_req(url,parse=True):
        '''Performs an HTTP request.
        
        Parameters
        ----------
        url: str
            URL for the request.
        parse: bool (optional)
            Whether to parse the result as JSON. Default: True.
        
        Returns
        -------
        result:
            Either a dictionary (if the request returns JSON and parse is True)
            or the raw result of the request.
        '''
        
        http = urllib3.PoolManager()
        r = http.request('GET', url)
        if parse:
            result = json.loads(r.data)
        else:
            result = r.data

        return result

    @classmethod
    def get_genes_in_class(cls,gene_class):
        '''Get genes (e.g. rab-1, rab-2, ...) in gene class (e.g. rab)
        
        Parameters
        ----------
        gene_class: str
            Gene class (e.g. rab is the class of the genes rab-1, rab-2, ...)
            
        Returns
        -------
        genes: list of str
            Genes in the gene class.
        '''
        
        url = "http://rest.wormbase.org/rest/field/gene_class/"+\
               gene_class+"/current_genes"
        result = cls.http_req(url)
        genes = result["current_genes"]["data"]["Caenorhabditis elegans"]
        
        return genes 
    
    @classmethod 
    def get_transcripts_ids(cls,gene_wbid):
        '''Get transcript isoform IDs for a specific gene. Does not check status
        (e.g. confirmed by cDNA).
        
        Parameters
        ----------
        gene_wbid: str or int
            WormBase ID of the gene. If string, the WBGene00000001 format is
            expected. 
            
        Returns
        -------
        t_ids: list of str
            List of transcripts for that gene. 
        
        '''
        
        if type(gene_wbid)==int:
            gene_wbid = cls.wbid_to_wbstr(gene_wbid)
        
        url = "http://rest.wormbase.org/rest/widget/gene/"+\
               gene_wbid+"/sequences"
        result = cls.http_req(url)
        if result["fields"]["name"]["data"]["taxonomy"] != "c_elegans":
            return None
        result = result["fields"]["gene_models"]["data"]["table"]
        t_ids = []
        for r in result:
            if "Coding transcript" in r["type"]: 
                for m in range(len(r["model"])):
                    t_id = r["model"][m]["id"]
                    t_ids.append(t_id)
                        
        return t_ids
    
    @classmethod
    def get_sequences_from_transcript_id(cls, transcript_wbid):
        url = "http://rest.wormbase.org/rest/widget/transcript/"+\
               transcript_wbid+"/sequences"
        r = cls.http_req(url)
        
        return r
        
    @classmethod
    def get_unspliced_sequence_from_transcript_id(cls, transcript_wbid):
        r = cls.get_sequences_from_transcript_id(transcript_wbid)["fields"]
        
        if r["strand"]["data"]=="-": strand = "negative_strand"
        else: strand = "positive_strand"
        
        r = r["unspliced_sequence_context_with_padding"]["data"]
        sequence = r[strand]["sequence"]
        features = r[strand]["features"]
        
        return {"sequence": sequence.lower(), "features": features}
        
