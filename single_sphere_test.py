import numpy as np
import owntools
import c2raytools as c2t
import matplotlib.pyplot as plt
import noise_kanan
import SLIC
import stitch_size_estimate

#z = 7.263

model_xf = np.zeros((250,250,250))
owntools.put_sphere(model_xf, [125,125,125], 50)

model_labels = SLIC.slic_cube(model_xf, n_segments=1000)
model_bin    = stitch_size_estimate.stitch_maximumdeviation(model_xf, model_labels, bins=50, binary=True)


