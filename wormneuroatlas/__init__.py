'''
Neural signal propagation atlas [1], genome [2], and single-cell 
transcriptome [3], neuropeptide/GPCR deorphanization [4], anatomical 
connectome [5,6], and monoaminergic connectome [7] all in one place. Worm Neuro 
Atlas allows to build a basic version of the neuropeptidergic connectome [8] 
([8] contains more detailed analysis). Please cite the code if you use it, along
with the papers containing the datasets you use.


1. Randi et al., arXiv 2022 https://arxiv.org/abs/2208.04790
2. WormBase, wormbase.org
3. Taylor et al., Cell 2021
4. Beets et al., bioRxiv 2022 https://www.biorxiv.org/content/10.1101/2022.10.30.514428v1
5. White et al., Phil. Trans. R. Soc 1986
6. Witvliet et al., Nature 2021
7. Bentley et al., PLOS Comp. Bio. 2016
8. Ripoll-Sanchez et al., bioRxiv 2022 https://www.biorxiv.org/content/10.1101/2022.10.30.514396v2
'''

__all__ = ["NeuroAtlas","WormBase","Cengen","PeptideGPCR","ExponentialConvolution_min"]

from .NeuroAtlas import NeuroAtlas
try:
    from .WormBase import WormBase
except:
    print("WormBase not available. Likely because libcurl is missing (see "+\
          "error below for details). To install libcurl, see the README)")
from .Cengen import Cengen
from .PeptideGPCR import PeptideGPCR
from .ExponentialConvolution_min import ExponentialConvolution_min
