import wormneuroatlas as wa

watlas = wa.NeuroAtlas()
                                    
gexp = watlas.get_gene_expression(gene_names=["npr-1"],th=4)
ASGL = watlas.ids_to_ais("ASGL")
print(gexp[ASGL])

gexp = watlas.cengen.get_expression(gene_wbids=["WBGene00006811"],neuron_ids=["ADA","ADE","ADF","ADL","AFD"],th=2)
print(gexp)
