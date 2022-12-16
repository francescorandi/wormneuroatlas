## Using the preprocessed, thresholded data from Alexis Weinreb's repository
## https://github.com/AlexWeinreb/cengenDataSC.

#install.packages("BiocManager")
#BiocManager::install("rhdf5")
library(rhdf5)

datafolder <- "/home/francesco/data/cengen/cengenDataSC-master/data-raw/"
dst <- "/home/francesco/dev/wormneuroatlas/wormneuroatlas/data/"

load(paste(datafolder,"cengen_nCells.rda",sep=""))
load(paste(datafolder,"cengen_proportion.rda",sep=""))
load(paste(datafolder,"cengen_sc_1.rda",sep=""))
load(paste(datafolder,"cengen_sc_2.rda",sep=""))
load(paste(datafolder,"cengen_sc_3.rda",sep=""))
load(paste(datafolder,"cengen_sc_4.rda",sep=""))
load(paste(datafolder,"cengen_TPM.rda",sep=""))
load(paste(datafolder,"cengen_UMI.rda",sep=""))
dim(cengen_sc_2)
dim(cengen_TPM)

h5fname <- paste(dst,"cengen.h5",sep="")
h5createFile(h5fname)

h5write(dimnames(cengen_TPM)[[1]],h5fname,"gene_wbids")
h5write(dimnames(cengen_TPM)[[2]],h5fname,"neuron_ids")
h5write(cengen_TPM,h5fname,"cengen_TPM")
h5write(cengen_nCells,h5fname,"cengen_nCells")
h5write(cengen_proportion,h5fname,"cengen_proportion")
h5write(cengen_sc_1,h5fname,"cengen_sc_1")
h5write(cengen_sc_2,h5fname,"cengen_sc_2")
h5write(cengen_sc_3,h5fname,"cengen_sc_3")
h5write(cengen_sc_4,h5fname,"cengen_sc_4")
h5write(dimnames(cengen_sc_1)[[2]],h5fname,"cengen_sc_1_neuron_ids")
h5write(dimnames(cengen_sc_1)[[1]],h5fname,"cengen_sc_1_gene_wbids")
h5write(dimnames(cengen_sc_2)[[2]],h5fname,"cengen_sc_2_neuron_ids")
h5write(dimnames(cengen_sc_2)[[1]],h5fname,"cengen_sc_2_gene_wbids")
h5write(dimnames(cengen_sc_3)[[2]],h5fname,"cengen_sc_3_neuron_ids")
h5write(dimnames(cengen_sc_3)[[1]],h5fname,"cengen_sc_3_gene_wbids")
h5write(dimnames(cengen_sc_4)[[2]],h5fname,"cengen_sc_4_neuron_ids")
h5write(dimnames(cengen_sc_4)[[1]],h5fname,"cengen_sc_4_gene_wbids")
h5write(cengen_UMI,h5fname,"cengen_UMI")
