import numpy as np
from hcipy.optics import Wavefront

def PSI(ref_, img_, I_sum, phi_sum, phi_2, phi_I, filter_, grid, prop_, i_):#, ncpa_filter_, normalisation_):
	#
	phi_I += ref_ * img_
	phi_2 += np.abs(ref_)**2
	phi_sum += ref_
	I_sum += img_

	psi_estimate_ = (phi_I - phi_sum * I_sum / (i_+1)) / (phi_2 - np.abs(phi_sum / (i_+1))**2)

	wf = Wavefront(psi_estimate_ * filter_)

	wf.electric_field *= np.exp(-2j * grid.as_('polar').theta)	# This line solves the phase wrapping
	pup = prop_.backward(wf)
	ncpa_estimate_ = pup.electric_field.imag

	return ncpa_estimate_, I_sum, phi_sum, phi_2, phi_I, psi_estimate_
