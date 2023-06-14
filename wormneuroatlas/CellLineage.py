import wormneuroatlas as wa
import numpy as np

class CellLineage:
    
    fname = "cell_lineage.txt"
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
        
        
    def load_from_file(self):
        f = open(self.module_folder+self.fname,"r")
        lines = f.readlines()
        f.close()
        
        self.cells = []
        self.lineages = []
        for l in lines:
            s = l.split("\t")
            cell = s[0]
            if "/" in cell: continue
            self.cells.append(cell)
            self.lineages.append(s[1])
            
    def get_lineage(self,cell):
        if cell in self.cells:
            return self.lineages[self.cells.index(cell)]
        else:
            return None
            
    def lineage_similarity(self, cell1, cell2, mode="end", normalize=True,
                           return_all=False):
        lin1 = self.get_lineage(cell1)
        lin2 = self.get_lineage(cell2)
        if lin1 is None or lin2 is None:
            return None
            
        spl_char1 = " " if "." not in lin1 else "."
        spl_char2 = " " if "." not in lin2 else "."
        founder1, lin1 = lin1.split(spl_char1)
        founder2, lin2 = lin2.split(spl_char2)
        
        same_founder = founder1==founder2
        
        lin1 = [s for s in lin1]
        lin2 = [s for s in lin2]
        
        n_lin1 = len(lin1)
        n_lin2 = len(lin2)
        
        matches = []
        ls = 0
        
        if mode=="all":
            if n_lin1!=n_lin2:
                raise ValueError("Cell lineages cannot be compared with"+\
                                 "mode=\"all\" if they don't have the same"+\
                                 "length.")
            
            for i in range(n_lin1):
                if lin1[i]==lin2[i]:
                    matches.append(True)
                    ls += 1
            
            if normalize:
                ls = ls/n_lin1
        
        elif mode=="end":
            n_lin12 = min(n_lin1,n_lin2)
            for i in range(n_lin12):
                if lin1[-i-1]==lin2[-i-1]: ls += 1
                else: break
            
            if normalize: ls = ls/n_lin12
            
        elif mode=="symmetric":
            if n_lin1!=n_lin2:
                raise ValueError("Cell lineages cannot be compared with"+\
                                 "mode=\"symmetric\" if they don't have the "+\
                                 "same length.")
            for i in range(n_lin1):
                if (lin1[i]=="a" and lin2[i]=="p") or\
                   (lin1[i]=="p" and lin2[i]=="a") or\
                   (lin1[i]=="l" and lin2[i]=="r") or\
                   (lin1[i]=="r" and lin2[i]=="l"):
                    matches.append(True)
                    ls += 1
                
            if normalize: ls = ls/n_lin1
            
        elif mode=="classify":
            if n_lin1!=n_lin2:
                raise ValueError("Cell lineages cannot be compared with"+\
                                 "mode=\"classify\" if they don't have the "+\
                                 "same length.")
            if n_lin1<2:
                ci = 0
            else:
                ci = 1
            
            if lin1[ci] in ["l","r"] and lin2[ci] in ["l","r"] and\
               all([lin1[ci+1+i]==lin2[ci+1+i] for i in range(n_lin1-ci-1)]):
                ls = 1
            else:
                ls = 0
            
        if return_all:
            return ls, same_founder
        else:
            return ls
        
    def lineage_similarity_old(self, cell1, cell2, c0="end", normalize=True):
        lin1 = self.get_lineage(cell1)
        lin2 = self.get_lineage(cell2)
        if lin1 is None or lin2 is None:
            return None
        
        if c0=="end":
            n_lin1 = len(lin1)
            n_lin2 = len(lin2)
            n_lin12 = min(n_lin1,n_lin2)
            n_match = 0
            for i in range(n_lin12):
                if lin1[-i]==lin2[-i]: n_match += 1
                else: break
            
            if normalize: n_match /= n_lin12
        
        return n_match
