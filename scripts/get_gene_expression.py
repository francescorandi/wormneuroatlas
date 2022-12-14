import wormneuroatlas as wa

watlas = wa.NeuroAtlas()

gexp = watlas.cengen.get_expression(gene_names=["avr-14","gur-3"],
                                    neuron_ids=["AVA","AVD"],
                                    threshold=1)

print(gexp)
