import numpy as np
import cupy as cp
import cupyx.scipy.ndimage as cndi
v=cp.asarray(np.random.random((48,256,256)).astype(np.float32))
print("uniform", float(cndi.uniform_filter(v,size=7).sum()))
print("median", float(cndi.median_filter(v,size=3).sum()))
print("zoom", cndi.zoom(v,(1.2,1,1),order=1).shape)
print("OK")
