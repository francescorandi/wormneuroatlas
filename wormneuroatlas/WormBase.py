import pycurl, certifi, json, numpy as np
from io import BytesIO

class WormBase:
    
    gene_ids_fname = "c_elegans.PRJNA13758.WS286.geneIDs.txt"
    
    def __init__(self):
        self.load_gene_ids()

    def load_gene_ids(self):
        '''Load gene IDs, names, other codes, and types from file. See 
        WormBase.gene_ids_fname for which version of the database is being 
        used.
        '''
        f = open(self.gene_ids_fname,"r")
        lines = f.readlines()
        f.close()
        self.g_ids = []
        self.g_names = []
        self.g_codes = []
        self.g_types = []
        for l in lines:
            la = l.split(",")
            if la[4] != "Live": continue
            self.g_ids.append(la[1])
            self.g_names.append(la[2])
            self.g_codes.append(la[3])
            self.g_types.append(la[5])
    
    
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

    @staticmethod
    def get_genes_in_class(gene_class):
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
        result = curl_req(url)
        genes = result["current_genes"]["data"]["Caenorhabditis elegans"]
        
        return genes 
    
    @staticmethod 
    def get_transcripts_ids(gene_id):
        '''Get transcript isoform IDs for a specific gene. Does not check status
        (e.g. confirmed by cDNA).
        
        Parameters
        ----------
        gene_id: str
            WormBase ID of the gene (e.g. WBGene00000001)
            
        Returns
        -------
        t_ids: list of str
            List of transcripts for that gene. 
        
        '''
        url = "http://rest.wormbase.org/rest/widget/gene/"+\
               gene_id+"/sequences"
        result = curl_req(url)
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
    
    @staticmethod
    def get_n_introns(transcript_id, starting_feature=None):
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
               transcript_id+"/sequences"
        r = curl_req(url)["fields"]
        if r["strand"]["data"]=="-": strand = "negative_strand"
        else: strand = "positive_strand"
        
        features = r["unspliced_sequence_context"]["data"][strand]["features"]
        feature_found = False
        if starting_feature is None:
            feature_found = True
            starting_feature = ""
        n_introns = 0
        for feature in features:
            if feature["type"] == :
                feature_found = True
            if feature_found and feature["type"]=="intron":
                n_introns += 1
        
        return n_introns
