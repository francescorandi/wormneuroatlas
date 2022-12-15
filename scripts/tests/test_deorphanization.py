import wormneuroatlas as wa

pg = wa.PeptideGPCR()

flp_1_gpcrs = pg.get_gpcrs_binding_to(["flp-1","flp-2"])
dmsr_1_peptides = pg.get_peptides_binding_to("dmsr-1")
print(flp_1_gpcrs)
print(dmsr_1_peptides) 
