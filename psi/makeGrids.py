import numpy as np
from scipy.signal import convolve2d
from .psi_utils import gauss_2Dalt
from hcipy.aperture import circular_aperture, make_obstructed_circular_aperture
from hcipy.coronagraphy import VortexCoronagraph
from hcipy.math_util import inverse_tikhonov
from hcipy.mode_basis import ModeBasis, make_zernike_basis
from hcipy.optics import Apodizer, make_gaussian_influence_functions, OpticalSystem

def makeFilters(grid, type="back_prop", sigma=0.05, cent_obs=0.27, outer_vin=1.0):
	"""Packages required:
	Circular and make_obstructed_circular_aperture
	gauss_2Dalt
	convolve2d
	"""
	if type == "back_prop":
		filter_ = circular_aperture(15)(grid)	# 15 is an arbitrary number set by Emiel.
	elif type == "ncpa":
		filter_ = make_obstructed_circular_aperture(1*0.7, 1.8*cent_obs)(grid)
	elif type == "reset":
		filter_ = make_obstructed_circular_aperture(1*outer_vin, 1.8*cent_obs)(grid)
	else:
		print("Type of filter not given!!")
		quit()
	sigma_ = int(sigma*np.sqrt(filter_.shape))
	if sigma_ != 0:
		kernal_ = gauss_2Dalt(size_x=int(np.sqrt(filter_.shape)), sigma_x=sigma_)
		filter_.shape = (int(np.sqrt(filter_.shape)),int(np.sqrt(filter_.shape)))
		filter_ = convolve2d(filter_, kernal_, mode='same')
		filter_ = filter_.ravel()
	return filter_

def makeMatrices(grid_, ao_acts_, aperture_, reconstruction_normalisation_):
	ao_modes_ = make_gaussian_influence_functions(grid_, ao_acts_, 1.2 / ao_acts_)							# Create an object containing all the available DM pistons, 1.0 to
	ao_modes_ = ModeBasis([mode * aperture_ for mode in ao_modes_])
	transformation_matrix_ = ao_modes_.transformation_matrix
	reconstruction_matrix_ = inverse_tikhonov(transformation_matrix_, reconstruction_normalisation_)
	return transformation_matrix_, reconstruction_matrix_

def makeOpticalSystem(grid_):
	lyot_stop_ = make_obstructed_circular_aperture(0.98, 0.3)(grid_)
	coro_      = OpticalSystem([VortexCoronagraph(grid_, 2), Apodizer(lyot_stop_)])
	return coro_

def makeZerns(n_zerns, grid, start):
	zernike_modes_ = make_zernike_basis(n_zerns, 1, grid, start)
	reconstructor_zernike_ = inverse_tikhonov(zernike_modes_.transformation_matrix, 1e-3)
	modal_means_ = np.ones(n_zerns)
	return zernike_modes_, reconstructor_zernike_, modal_means_
