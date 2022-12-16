import wormneuroatlas as wa

pg = wa.PeptideGPCR()

flp_1_gpcrs = pg.get_gpcrs_binding_to(["pdf-1","flp-1"])
dmsr_1_peptides = pg.get_peptides_binding_to("ntr-1")
print(flp_1_gpcrs)
print(dmsr_1_peptides) 
