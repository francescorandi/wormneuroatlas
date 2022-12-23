import wormneuroatlas as wa

pg = wa.PeptideGPCR()

flp_1_gpcrs = pg.get_gpcrs_binding_to(["pdf-1","flp-1"])
frpr_8_peptides = pg.get_peptides_binding_to("frpr-8")

print(flp_1_gpcrs)
print(frpr_8_peptides) 
