import numpy as np, h5py, datetime
import wormneuroatlas as wa

atlas = wa.NeuroAtlas()

dst = "/home/francesco/dev/wormneuroatlas/wormneuroatlas/data/aconnectome_default.h5"
f = h5py.File(dst,"w")

time_compiled = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
f.attrs["time_compiled"] = np.string_(time_compiled)
f["chem"] = atlas.aconn_chem
f["gap"] = atlas.aconn_gap

f.close()
