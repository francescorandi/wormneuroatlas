# Worm Neuro Atlas
Neural signal propagation atlas [1], genome [2], and single-cell transcriptome [3], neuropeptide/GPCR deorphanization [4], anatomical connectome [5,6], and monoaminergic connectome [7] all in one place. Worm Neuro Atlas allows to build a basic version of the neuropeptidergic connectome [8] ([8] contains more detailed analysis). Please cite the code if you use it, along with the papers containing the datasets you use.

**Read the full documentation here** https://francescorandi.github.io/wormneuroatlas/

**Run a demo of Worm Neuro Atlas in a Jupyter notebook on Colab here** https://colab.research.google.com/drive/1j59Mv_PSaD4Nj2ITNDQWVxWfrl6ssqXy?usp=sharing

## Requirements
Interfacing with Worm Base requires the curl library. If the library is not present, the Worm Neuro Atlas will work but the Worm Base functionalities will not be available. On Windows, I was able to install and use `wormneuroatlas` without the need to install curl separately. If you encounter problems on Windows, let me know. On Debian/Ubuntu, libcurl can be installed via
```
sudo apt install libcurl4-openssl-dev libssl-dev
```

## Installation
`wormenuroatlas` can be installed via pip 
```
pip install wormneuroatlas
```
or by cloning this repository and pip install the module, e.g. via
```
git clone git@github.com:francescorandi/wormneuroatlas.git
cd wormneuroatlas
python -m pip install .
```


## NeuroAtlas class
NeuroAtlas is the main class that aggregates all the datasets, and directly handles the Signal Propagation Atlas, the anatomical connectome, and the monoaminergic connectome [7]. 

You can access the wild-type and unc-31 signal propagation atlas, for example, via
```
NeuroAtlas.get_signal_propagation_atlas(strain="wt") # or strain="unc31"
```

NeuroAtlas will use the other classes of the Python module (WormFunctionalConnectivity, WormBase, Cengen, PeptideGPCR, ...) to access data from other datasets like neural signal propagation, genome, single-cell transcriptome and aggregate those datasets. 
For example, you can get gene expression levels of genes flp-1 and aqp-1 in neurons AVAL and AVDR via
```
NeuroAtlas.get_gene_expression(gene_names=["flp-1","aqp-1"], neuron_names=["AVAL","AVDR"]).
```

See __Examples__ for more information on how to access the various datasets.

NeuroAtlas also manages the conversions between different conventions for neural IDs. NeuroAtlas can be instantiated to maintain the "exact" neural identities, or to merge neurons into classes (i.e. to approximate neuron identities): `merge_bilateral=True` will merge results for, e.g., AVAL and AVAR into the class AVA_, `merge_dorsoventral=True` will merge RMED and RMEV into RME_, while `merge_numbered=True` will merge VB3, VB4, ... into VB. These options can be combined to merge, for example, SMBVL, SMBVR, SMBDL, and SMBDR into SMB__ with `merge_bilateral=True, merge_dorsoventral=True`. 

## Cengen class
The Cengen class interfaces with single-cell RNASeq database from the CeNGEN project. Its main function is Cengen.get_gene_expression(), which is called by `NeuroAtlas.get_gene_expression()` after converting neuron IDs to CeNGEN-style IDs.

## PeptideGPCR
The PeptideGPCR class provides an interface to the neuropeptide/GPCR deorphanization in [4]. The two main functions are
```
PeptideGPCR.get_gpcrs_binding_to(peptides)
```
and
```
PeptideGPCR.get_peptides_binding_to(gpcrs)
```
which return the GPCRs binding to given peptides and the peptides binding to given GPCRs, respectively.

## WormBase class
The WormBase class uses the REST API provided by wormbase.org. WormBase currently has methods to retrieve lists of transcripts for given genes, as well as functions to convert WormBase-style gene IDs to gene names, etc.

## Examples
* plot_signal_propagation.py shows you how to access the signal propagation data,
* get_gene_expression.py shows you how to get gene expression data from CeNGEN,
* test_deorphanization.py shows you how to get the peptides binding to given GPCRs and vice versa,
* make_peptidergic_connectome.py shows you how to combine these functions to compile the neuropeptidergic connectome [8] using gene expression and neuropeptide/GPCR deorphanization.

## References
1. Randi et al., arXiv 2022 https://arxiv.org/abs/2208.04790
2. WormBase, wormbase.org
3. Taylor et al., Cell 2021
4. Beets et al., bioRxiv 2022 https://www.biorxiv.org/content/10.1101/2022.10.30.514428v1
5. White et al., Phil. Trans. R. Soc 1986
6. Witvliet et al., Nature 2021
7. Bentley et al., PLOS Comp. Bio. 2016
8. Ripoll-Sanchez et al., bioRxiv 2022 https://www.biorxiv.org/content/10.1101/2022.10.30.514396v2
