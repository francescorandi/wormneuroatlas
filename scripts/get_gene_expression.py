import wormneuroatlas as wa

watlas = wa.NeuroAtlas()

gexp = watlas.cengen.get_expression(gene_names=["avr-14","gur-3"],
                                    neuron_ids=["AVA","AVD","BAG"],
                                    threshold=1)
                                    
gexp2 = watlas.get_gene_expression(gene_names=["avr-14","gur-3"],
                                    neuron_ids=["AVAL","AVDR","BAGL"],
                                    threshold=1)

print(gexp)
print(gexp2)
