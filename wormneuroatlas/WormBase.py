import pycurl, certifi, json, numpy as np
from io import BytesIO
import wormneuroatlas as wa

class WormBase:
    
    gene_wbids_fname = "c_elegans.PRJNA13758.WS286.geneIDs.txt"
    module_folder = "/".join(wa.__file__.split("/")[:-1])+"/data/"
    '''Folder of the pumpprobe module'''
    
    def __init__(self, gene_subset = None):
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
    def curl_req(url):
        '''General cURL request to Wormbase. Look up the URLs to fetch specific
        data here: http://rest.wormbase.org/index.html .
        
        Parameters
        ----------
        url: str
            URL of the request.
        
        Returns
        -------
        result: dictionary
            JSON-parsed result of the cURL request.
        
        '''
        # Creating a buffer as the cURL is not allocating a buffer for the network response
        buffer = BytesIO()
        c = pycurl.Curl()
        #initializing the request URL
        c.setopt(c.URL, url)
        #setting options for cURL transfer  
        c.setopt(c.WRITEDATA, buffer)
        #setting the file name holding the certificates
        c.setopt(c.CAINFO, certifi.where())
        # perform file transfer
        c.perform()
        #Ending the session and freeing the resources
        c.close()

        #retrieve the content BytesIO
        result = json.loads(buffer.getvalue().decode('iso-8859-1'))
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
        result = cls.curl_req(url)
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
        result = cls.curl_req(url)
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
    def get_n_introns(cls,transcript_wbid, starting_feature=None):
        '''Get the number of introns for a given transcript isoform. If 
        starting_feature is specified, counts introns only starting from that
        feature in the gene.
        
        Parameters
        ----------
        transcript_id: str
            ID of the transcript.
        starting_feature: str (optional)
            Feature where to start counting. Default: None.
            
        Returns
        -------
        n_introns: int
            Number of introns.        
        '''
    
        url = "http://rest.wormbase.org/rest/widget/transcript/"+\
               transcript_wbid+"/sequences"
        r = cls.curl_req(url)["fields"]
        if r["strand"]["data"]=="-": strand = "negative_strand"
        else: strand = "positive_strand"
        
        features = r["unspliced_sequence_context"]["data"][strand]["features"]
        feature_found = False
        if starting_feature is None:
            feature_found = True
            starting_feature = ""
        n_introns = 0
        for feature in features:
            if feature["type"] == starting_feature:
                feature_found = True
            if feature_found and feature["type"]=="intron":
                n_introns += 1
        
        return n_introns
