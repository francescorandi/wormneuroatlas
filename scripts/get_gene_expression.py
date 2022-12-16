import wormneuroatlas as wa

watlas = wa.NeuroAtlas()
                                    
#gexp = watlas.get_gene_expression(gene_names=["pdf-1"],threshold=4)
#gexp = watlas.get_gene_expression(neuron_ids=["AQR"],gene_names=["pdf-1"],threshold=4)
#print(gexp)
gexp = watlas.get_gene_expression(gene_names=["npr-1"],th=3)
ASGL = watlas.ids_to_ais("ASGL")
print(gexp[ASGL])
gexp = watlas.cengen.get_expression(gene_names=["npr-1"],neuron_ids=["ASG"],th=3)
print(gexp)

gexp = watlas.cengen.get_expression(gene_wbids=["WBGene00006811"],neuron_ids=["ADA","ADE","ADF","ADL","AFD"],th=2)
print(gexp)
