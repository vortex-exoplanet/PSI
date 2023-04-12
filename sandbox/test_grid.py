'''
testing grids with and without physical dimensions
'''

import hcipy
import numpy as np
import matplotlib.pyplot as plt

ngrid =128
det_size=15
det_res=4
diam = 40

pupilGrid = hcipy.make_pupil_grid(ngrid, diameter=diam)
focalGrid = hcipy.make_focal_grid(det_res,
                                  det_size,
                                #   spatial_resolution= diam,
                                  pupil_diameter=diam,
                                  reference_wavelength=1,
                                  focal_length=diam)

prop = hcipy.FraunhoferPropagator(pupilGrid, focalGrid, focal_length=diam)

optical_model = hcipy.OpticalSystem([prop])

aperture = hcipy.make_circular_aperture(diam)(pupilGrid)

phase = hcipy.Field(0, pupilGrid)
wf = hcipy.Wavefront(np.exp(1j * phase) * aperture)

image = optical_model(wf).power

plt.figure()
hcipy.imshow_field(image)


import psi.psi_utils as psi_utils
filter_fp = psi_utils.makeFilters(focalGrid,
                                  "back_prop",
                                 sigma=0.1,
                                lD = 6 * 2)

sh = int(np.sqrt(len(filter_fp)))
plt.figure()
plt.imshow(filter_fp.reshape(sh,sh))

# psi_sensor.filter_fp = psi_utils.makeFilters(psi_sensor.inst.focalGrid,
#                                   "back_prop",
#                                  sigma=0.05,
#                                 lD = 2 * 2)