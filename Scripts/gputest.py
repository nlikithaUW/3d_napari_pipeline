"""Smoke-test that cupy can actually run the ndimage kernels this pipeline uses.

Run this first when the GPU path misbehaves — it isolates "is CUDA/cupy working
at all" from anything in `pipeline`. A working install prints four lines and
"OK"; a broken one raises here rather than silently degrading to CPU inside
`pipeline.gpu`.

    python Scripts/gputest.py
"""

import numpy as np
import cupy as cp
import cupyx.scipy.ndimage as cndi
v=cp.asarray(np.random.random((48,256,256)).astype(np.float32))
print("uniform", float(cndi.uniform_filter(v,size=7).sum()))
print("median", float(cndi.median_filter(v,size=3).sum()))
print("zoom", cndi.zoom(v,(1.2,1,1),order=1).shape)
print("OK")
