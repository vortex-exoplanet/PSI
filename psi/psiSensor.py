

import time
import os
import sys
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
import getpass
import datetime
import shutil

import hcipy
import importlib

# sys.path.append('/Users/orban/Projects/METIS/4.PSI/psi_github/')
import psi.psi_utils as psi_utils
from .configParser import loadConfiguration
from .instruments import  CompassSimInstrument, DemoCompassSimInstrument
from .helperFunctions import LazyLogger, timeit, build_directory_name, copy_cfgFileToDir

from astropy.visualization import imshow_norm,\
    SqrtStretch, MinMaxInterval, PercentileInterval, \
    LinearStretch, SinhStretch, LogStretch


# from config.config_metis_compass import conf

# sys.path.append('/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/')
print(hcipy.__file__)
print(hcipy.__version__)
# assert 1==0
# from hcipy import *
# import hcipy

class PsiSensor():
	'''
	 Phase Sorting Interferometry wavefront sensor

	 Parameters
	 ----------
	 config_file : str
	 	filename of the Python config file
	 logger : object
	 	logger object. Default is ``LazyLogger``
	'''
	def __init__(self, config_file, logger=LazyLogger('PSI')):
		self.logger=logger
		self.logger.info('Loading and checking configuration')
		self._config_file = config_file
		self.cfg = loadConfiguration(config_file)


	def setup(self):
		'''
			Setup the PSI wavefront sensor based on the configuration
		'''
		# Build instrument object 'inst'
		self.logger.info('Initialize the instrument object & building the optical model')
		# self.inst = getattr(instruments,
		# 					self.cfg.params.instrument)(self.cfg.params)
		self.inst = eval(self.cfg.params.instrument)(self.cfg.params)
		importlib.import_module
		self.inst.build_optical_model()

		# # Build focal plane filter for PSI
		self.filter_fp = psi_utils.makeFilters(self.inst.focalGrid,
		                                     "back_prop",
		                                     sigma=self.cfg.params.psi_filt_sigma,
		                                     lD = self.cfg.params.psi_filt_radius * 2)



		# Build basis for modal projection
		# TODO move to configParser compute ???
		# TODO orthogonalization should be done on the specific aperture
		if self.cfg.params.psi_correction_mode is not 'all':
			self.logger.info('Building modal basis for projectiong/filtering of the NCPA map')
			diam = 1
			if self.cfg.params.psi_correction_mode == 'zern':
			    self.M2C = hcipy.make_zernike_basis(self.cfg.params.psi_nb_modes, diam,
										 	   self.inst.pupilGrid,
										 	   self.cfg.params.psi_start_mode_idx)#.orthogonalized
			if self.cfg.params.psi_correction_mode == 'dh':
				self.logger.warn('Warning psi_start_mode_idx is ignored')
				self.M2C = hcipy.make_disk_harmonic_basis(self.inst.pupilGrid,
													 self.cfg.params.psi_nb_modes,
													 diam)

			self.C2M = hcipy.inverse_tikhonov(self.M2C.transformation_matrix, 1e-3)

		self._ncpa_modes_integrated = 0 # np.zeros(self.cfg.params.psi_nb_modes)

		# Diffraction component to be removed from speckle field in pupil
		self._diffraction_component= self.inst.aperture

		# Holders for PSI-related cubes
		N = int(np.round(1 / (self.cfg.params.dit * self.cfg.params.psi_framerate)))

		self._speckle_field_t_psi = np.zeros((N, self.inst.focalGrid.shape[0] *\
											  self.inst.focalGrid.shape[1]),
									         dtype=np.complex_)
		self._image_t_psi= np.zeros((N, self.inst.focalGrid.shape[0] *\
					           			self.inst.focalGrid.shape[1]))


		# Photometry
		self.nbOfPhotons = self.inst.getNumberOfPhotons()

		# -- [WIP] Mask for NCPA correction
		# if cfg.params.inst_mode == 'CVC':
		# 	mask = inst.lyot_stop_mask
		# else:
		# 	mask = inst.aperture
		self.ncpa_mask = self.inst.aperture



		# -- [WIP] Scaling of ncpa correction
		self.ncpa_scaling =  1 #1.e6  # NB: this is not the same as cfg.params.ncpa_scaling !!
		self._ncpa_correction_long_term = 0 # 24/06/2022 -- for WV integrating

		self._skip_limit = self.cfg.params.psi_skip_limit

		self.iter = 0 # iteration index of PSI
		self._loop_stats = []

		# -- Plotting & saving results
		# self.fig = plt.figure(figsize=(9, 3))
		if self.cfg.params.save_loop_statistics:
			self._directory = build_directory_name(self._config_file,
											 self.cfg.params.save_basedir)
			if not os.path.exists(self._directory):
				os.makedirs(self._directory)

			# copy config file to directory
			copy_cfgFileToDir(self._directory, self._config_file)

			if self.cfg.params.save_phase_screens:
				self._directory_phase = self._directory + 'residualNCPA/'
				os.mkdir(self._directory_phase)

			self.logger.info('Results will be stored in {0}'.format(self._directory))


	def _save_loop_stats(self):
		'''
			Saving loop statistics to file
		'''
		data = np.array(self._loop_stats)
		np.savetxt(os.path.join(self._directory, 'loopStats.csv'),
				  data,
				  header ='units are nm \n it \t wfe_all_f \t wfe_qs_f \t wfe_all \t wfe_qs',
				  fmt=['%i' , '%f', '%f', '%f', '%f'],
				  delimiter= '\t')

	def _store_phase_screens_to_file(self, i):
		'''
			Storing phase screens to fits file.
			Units of phase is nm

			TODO populate the header with useful information
		'''
		conv2nm = self.inst.wavelength / (2*np.pi) * 1e9

		ncpa_correction = self.inst.phase_ncpa_correction * conv2nm
		ncpa_injected = (self.inst.phase_ncpa + self.inst.phase_wv) * conv2nm

		if type(ncpa_correction) == hcipy.Field:
			ncpa_correction = np.copy(ncpa_correction.shaped)
		if type(ncpa_injected) == hcipy.Field:
			ncpa_injected = np.copy(ncpa_injected.shaped)

		ncpa_residual = ncpa_injected + ncpa_correction

		filename = 'residual_ncpa_%s.fits'%i
		full_name = self._directory_phase + filename
		fits.writeto(full_name, ncpa_residual)
		hdr = fits.getheader(full_name)
		hdr.set('EXTNAME', 'NCPA_IN')
		fits.append(full_name, ncpa_injected, hdr)
		hdr = fits.getheader(full_name)
		hdr.set('EXTNAME', 'NCPA_COR')
		fits.append(full_name, ncpa_correction, hdr)


	def _psiCalculation(self, speckle_fields_fp, images_fp, scale_factor=False):
		'''
		1. Calculate the 'subject beam' :math:`\Psi`, which is the electric field in the \
			focal plane corresponding to the NCPA we want to estimate. See eq. 21 in Codona et all. 2017:

		 .. math::
				\Psi = \dfrac{<\psi I >}{<|\psi^2|>}

		2. Then propagate backwards to the pupil plane.

		3. Based on the small aberration hypothesis, return the imaginary part\
		   as the NCPA phase map estimate

		N.B.: speckle_fields_fp & images_fp should have the same dimension and be in sync.

		Parameters
		----------
		speckle_fields_fp : array
			speckle field in the focal plane calculated based on the WFS telemetry
		images_fp : array
			science images

		'''
		# PSI calculation
		phi_I = np.sum(speckle_fields_fp * images_fp, axis=0)
		phi_2 = np.sum(np.abs(speckle_fields_fp)**2, axis=0)
		phi_sum = np.sum(speckle_fields_fp, axis=0)
		I_sum = np.sum(images_fp, axis=0)
		nbframes = images_fp.shape[0]

		# ff = np.sum(I_sum) / np.sum(phi_2)
		# phi_2 *= ff
		# phi_I *=np.sqrt(ff)
		# phi_sum *= np.sqrt(ff)
		#---------
		if scale_factor:
			# Trying to compute the scale factor
			var_I = np.var(images_fp, axis=0)
			var_s = np.var(speckle_fields_fp**2, axis=0)
			# s_mean = np.mean(np.abs(speckle_field_t_psi)**2, axis=0)

			phi_I_c = (phi_I - phi_sum * I_sum / nbframes)
			phi_2_c = (phi_2 - np.abs(phi_sum / nbframes)**2)

			I_mean = np.mean(images_fp, axis=0)
			g = (var_I / var_s - 2 * np.abs(phi_I)**2 / (phi_2 * var_s))**(1/4)


		else:
			g = 1
		#---------

		psi_estimate = (phi_I - phi_sum * I_sum / nbframes) / (g * (phi_2 - np.abs(phi_sum / nbframes)**2))
		# psi_estimate = (phi_I) / (g * (phi_2 ))


		self._psi_estimate = hcipy.Field(psi_estimate, self.inst.focalGrid)
		wf = hcipy.Wavefront(self._psi_estimate * self.filter_fp)

		# Propagate estimation back to the entrance pupil
		pup = self.inst.optical_model.backward(wf)
			## -- alternative for the SVC charge 2 --
			# wf.electric_field *= np.exp(-2j * inst.focalGrid.as_('polar').theta)
			# pup = inst._prop.backward(wf)

		self._estimated_wavefront = pup

		# Small phase hypothesis
		ncpa_estimate = pup.electric_field.imag

		return ncpa_estimate

	def _projectOnModalBasis(self, ncpa_estimate):
		'''
			Projection of NPCA phase map onto a finite set of modes.
			Projection matrices are defined in self.setup

			Parameters
			---------
			ncpa_estimate	: NCPA phase map

			Return
			------
			ncpa_estimate	: phase map filtered to a finite set of modes
			ncpa_modes		: modal coefficients vector

		'''
		# Project ncpa estimate on finite set of modes
		if self.cfg.params.inst_mode == 'CVC' or self.cfg.params.inst_mode == 'RAVC':
			proj_mask = self.inst.lyot_stop_mask
		else:
			proj_mask = self.inst.aperture
		ncpa_modes      = self.C2M.dot(ncpa_estimate * proj_mask)
		ncpa_estimate  = self.M2C.transformation_matrix.dot(ncpa_modes)
		return ncpa_estimate, ncpa_modes

	def _propagateSpeckleFields(self, wfs_telemetry_buffer):
		'''
		Compute a focal-plane complex speckle fields using WFS telemetry

		Parameters
		----------
			wfs_telemetry_buffer : cube of WFS telemetry phase map

		Returns
		-------
			speckle_fiels : corresponding complex speckle fields
		'''
		nf, nx, ny = wfs_telemetry_buffer.shape
		wfs_telemetry_buffer_1d = wfs_telemetry_buffer.reshape((nf, nx*ny))
		wfs_wavefront_hcipy = hcipy.Field(wfs_telemetry_buffer_1d, self.inst.pupilGrid)
		Efield = hcipy.Wavefront(self.inst.aperture * np.exp(1j * wfs_wavefront_hcipy) \
			- self._diffraction_component)
		Efield.total_power = self.nbOfPhotons  * nf

		# #----------
		# # [2022-06-22] revising the flux scaling
		# # flux-perfect and flux_speckle could be computed only once
		# # this does not work yet ...
		# Efield_perfect = hcipy.Wavefront(self._diffraction_component)
		# flux_perfect = np.copy(Efield_perfect.total_power)
		# Efield_speckle = hcipy.Wavefront(self.inst.aperture * 1j * wfs_wavefront_hcipy[0] )
		# flux_speckle = np.copy(Efield_speckle.total_power)
		#
		# Efield = hcipy.Wavefront(self.inst.aperture * np.exp(1j * wfs_wavefront_hcipy) \
		# 	- self._diffraction_component)
		# Efield.total_power = self.nbOfPhotons  * nf * (flux_speckle / flux_perfect)
		# # Efield = hcipy.Wavefront(Efield_telemetry.electric_field - Efield_perfect.electric_field)
		# self.logger.debug('Efield total power= {0} vs Efield_perfect = {1}, Efield_telemetry = {2} '.format(Efield.total_power,
		# 	Efield_perfect.total_power, Efield_speckle.total_power))
		# #--------


		# Efield_perfect = hcipy.Wavefront(self._diffraction_component)
		# Efield_perfect.total_power = self.nbOfPhotons / nf
		#
		speckle_fields = self.inst.optical_model(Efield)
		# speckle_fields.total_power = self.nbOfPhotons * nf
		# speckle_fields_perfect = self.inst.optical_model(Efield_perfect).electric_field
		#
		# speckle_fields = speckle_fields - speckle_fields_perfect
		# reproducing what is normally done in Wavefront.power to obtain the correct flux normalisation
		return speckle_fields.electric_field * np.sqrt(speckle_fields.grid.weights)

	@timeit
	def _fullPsiAlgorithm(self, wfs_telemetry_buffer, science_images_buffer):
		'''
			Complete PSI algorithm containing the following steps:
			1. synchronization of the WFS phase telemetry and the science images buffers
			2. Computation of the speckle fields based on the WFS telemetry.
				This corresponds to the "reference beam" :math:`\psi`
			3. PSI algebra using the :math:`\psi` buffer, and the :math:`I` (images) buffer
			4. Optional projection on a finite set of modes

			Parameters
			----------
			wfs_telemetry_buffer

			science_images_buffer

			Returns
			-------
			ncpa_estimate
		'''
		# Synchronize the WFS telemetry buffer and science image buffer
		# Telemetry_indexing is used for sync and slicing of the wfs telemetry buffer
		telemetry_indexing = self.inst.synchronizeBuffers(wfs_telemetry_buffer,
													      science_images_buffer)

		# PSI : compute speckle fields in the focal plane from the telemetry buffer
		# speckle_fields = self._propagateSpeckleFields(wfs_telemetry_buffer)

		# Slicing : one speckle field per science image
		for i in range(len(telemetry_indexing)):
			_s, _e = telemetry_indexing[i]
			# average_speckle_field = speckle_fields[_s:_e,:].sum(axis=0)
			speckle_fields = self._propagateSpeckleFields(wfs_telemetry_buffer[_s:_e])
			average_speckle_field = speckle_fields.mean(axis=0)

			self._speckle_field_t_psi[i] = average_speckle_field
			self._image_t_psi[i] = science_images_buffer[i].ravel()

		# TODO check if this conversion to hcipy.Field is needed
		self._speckle_field_t_psi = hcipy.Field(self._speckle_field_t_psi,
												self.inst.focalGrid)
		self._image_t_psi = hcipy.Field(self._image_t_psi,
										self.inst.focalGrid)

		# Raw PSI estimate
		ncpa_estimate = self._psiCalculation(self._speckle_field_t_psi,
										      self._image_t_psi)

		# Optional: project on modal basis
		if self.cfg.params.psi_correction_mode is not 'all':
			ncpa_estimate, ncpa_modes = self._projectOnModalBasis(ncpa_estimate)
			return ncpa_estimate, ncpa_modes
		else:
			return ncpa_estimate


	def findNcpaScaling(self, ncpa_estimate, rms_desired=None):
		'''
			Compute a NCPA scaling based on the input rms and
			the expected rms (provided in the config file).

			This scaling is later use to scale the PSI estimate before NCPA
			correction.

			Parameters
			---------
			ncpa_estimate

			Returns
			-------
			NCPA scaling
		'''
		conv2nm = self.inst.wavelength / (2 * np.pi) * 1e9
		rms_estimate = np.std(ncpa_estimate[self.ncpa_mask==1]) * conv2nm
		if rms_desired is None:
			rms_expected = self.cfg.params.ncpa_expected_rms
		else:
			rms_expected = rms_desired

		scaling = rms_expected / rms_estimate

		return scaling


	def next(self, display=True, check=False):
		'''
			Perform a complete iteration. This consists in:
			1. grab the WFS telemetry and the sciences image
			2. run the PSI algorithm
			3. set the NCPA correction
			4. (optional) check convergence
			5. (optional) show progress

			PARAMETERS
			-----------
			skip_limit	: float
				nm rms WFE above which the NCPA estimate will be skipped.
				For no skip_limit, use 'None'
		'''
		# Acquire telemetry buffers
		nbOfSeconds = 1/self.cfg.params.psi_framerate
		wfs_telemetry_buffer = self.inst.grabWfsTelemetry(nbOfSeconds)
		science_images_buffer = self.inst.grabScienceImages(nbOfSeconds)

		# Compute NCPA
		self._ncpa_estimate, self._ncpa_modes = self._fullPsiAlgorithm(wfs_telemetry_buffer,
										  	   science_images_buffer)

		if self.iter == 0:
			''' at first iteration, compute a NCPA scaling'''
			scaling = self.findNcpaScaling(self._ncpa_estimate)
			self.logger.info('New ncpa scaling is {0}'.format(scaling))
			self.ncpa_scaling = scaling

		# Arbitratry gain rule
		# For the first 5 iteration, this gives: [1.0, 0.5, 0.25, 0.125, 0.1]
		gain = np.max((0.8**self.iter, 0.45))
		# gain = 0.45 # 2022-06-24 --- dominated by water vapour


		ncpa_command = - gain * self._ncpa_estimate * self.ncpa_mask * self.ncpa_scaling

		self._ncpa_modes_integrated = self._ncpa_modes_integrated +\
		 	gain * self._ncpa_modes * self.ncpa_scaling

		if self._skip_limit is not None:
			ncpa_estimate_rms = np.sqrt(np.sum(self._ncpa_modes**2)) * \
			 	self.ncpa_scaling * self.inst.wavelength / 6.28 * 1e9
			# scaling = self.findNcpaScaling(ncpa_command, rms_desired=50)
			# print('Debug scaling : {0}'.format(scaling))
			if ncpa_estimate_rms > skip_limit :
				self.logger.warning('NCPA estimate too large ! Skipping !')
				ncpa_command= 0 * ncpa_command

		# Send correction
		self.inst.setNcpaCorrection(ncpa_command)

		self.iter += 1
		# ------------#
		# Inspect PSI convergence
		if check:
			self.checkPsiConvergence()

		# Metrics... - might only be valid in simualtion
		# self.evaluateSensorEstimate()

		# Display
		if display:
			I_avg = science_images_buffer.mean(0)
			self.show(I_avg,
					  self._ncpa_estimate * self.ncpa_scaling,
					  gain * self._ncpa_modes * self.ncpa_scaling)

	def loop(self):
		'''
			Run PSI for a number of iterations.
			At each iterations:
				1. run ``next()``
				2. evaluate the sensor estimate performance
				3. (optional) save fits file at every iteration
				4. (optional) save loop statistics at the end of the for-loop
		'''
		for i in range(self.cfg.params.psi_nb_iter):
			self.next()
			self.evaluateSensorEstimate()
			if self.cfg.params.save_phase_screens:
				self._store_phase_screens_to_file(self.iter)

		if self.cfg.params.save_loop_statistics:
			self._save_loop_stats()

	def show(self, I_avg, ncpa_estimate, ncpa_modes):
		'''
			Display the PSF and the NCPA correction

			TODO improve and add displays
		'''
		# self.fig.clf()
		# self.fig.gca()
		# plt.figure()
		plt.clf()
		gs = gridspec.GridSpec(2, 3)
		ax = plt.subplot(gs[0, 0])
		# hcipy.imshow_field(np.log10(I_sum / nbframes / I_sum.max()),
		#                    vmin=-4, vmax=-1.5)
		# hcipy.imshow_field(np.sqrt(I_sum / nbframes))
		# plt.imshow(np.sqrt(I_avg))
		im, norm = imshow_norm(I_avg, plt.gca(), origin='lower',
			interval=MinMaxInterval(),
			stretch=SqrtStretch())
		plt.axis('off')


		plt.title('Average SCI image')
		ax = plt.subplot(gs[0, 1])
		if np.size(self.inst.phase_ncpa_correction) == 1:
		    hcipy.imshow_field(np.zeros(256**2), self.inst.pupilGrid,
							   cmap='RdBu', mask=self.ncpa_mask)
		else:
		    # , vmin=-ncpa_max, vmax=ncpa_max)
		    hcipy.imshow_field(-self.inst.phase_ncpa_correction,
								self.inst.pupilGrid,
								cmap='RdBu', mask=self.ncpa_mask)
		plt.axis('off')
		plt.title('NCPA correction')

		ax = plt.subplot(gs[0, 2])
		hcipy.imshow_field(ncpa_estimate * self.ncpa_mask)
		plt.axis('off')
		plt.title('Last NCPA estimate')

		ax = plt.subplot(gs[1, :])
		if self.cfg.params.psi_correction_mode is not 'all':
			mm=np.arange(self.cfg.params.psi_start_mode_idx,
				self.cfg.params.psi_nb_modes + self.cfg.params.psi_start_mode_idx)
			plt.plot(mm, ncpa_modes, label='last NCPA correction')
			plt.plot(mm, self._ncpa_modes_integrated, c='k', ls='--', label='integrated')
			# plt.title('Last NCPA modes')
			plt.legend()
			plt.ylim((-0.1, 0.1))
			plt.xlabel('Mode index')
			plt.ylabel('rms [rad]')

		plt.draw()
		plt.pause(0.01)

	def checkPsiConvergence(self):
		'''
			TODO use self._speckle_field_t_psi and self._image_t_psi
				to compute the psiEstimate and see how it converge
		'''
		nbSteps= self._image_t_psi.shape[0]
		ncpa_estimates = np.zeros((nbSteps, self._ncpa_estimate.shape[0]))
		for i in range(1, nbSteps):
			tmp_estimate = self.ncpa_mask * self._psiCalculation(self._speckle_field_t_psi[:i,:],
											      self._image_t_psi[:i,:])

			if self.cfg.params.psi_correction_mode is not 'all':
				ncpa_estimate, ncpa_modes = self._projectOnModalBasis(tmp_estimate)
				ncpa_estimates[i,:] = ncpa_estimate
			else:
				ncpa_estimates[i,:] = ncpa_estimate

		return ncpa_estimates

	def evaluateSensorEstimate(self, verbose=True):
		'''
			Compute the rms errors made on quasi-static NCPA and on water vapour seeing.

			Only valid for a CompassSimInstrument and DemoCompassSimInstrument

			TODO make it generic to any instruments
		'''
		res_ncpa_qs = self.inst.phase_ncpa + self.inst.phase_ncpa_correction
		res_ncpa_all = self.inst.phase_ncpa + self.inst.phase_wv + \
			self.inst.phase_ncpa_correction
		if self.iter == 0:
			res_static_ncpa_qs = self.inst.phase_ncpa
		else:
			# tmp_avg = np.mean(self.inst.phase_ncpa_correction[self.inst.aperture>=0.5])
			self._ncpa_correction_long_term += self.inst.phase_ncpa_correction #- tmp_avg)
			# self._ncpa_correction_long_term /= self.iter

			res_static_ncpa_qs = self.inst.phase_ncpa + (self._ncpa_correction_long_term / self.iter)

		conv2nm = self.inst.wavelength / (2 * np.pi) * 1e9
		# rms_input_qs = np.std(self.inst.phase_ncpa[self.inst.aperture==1]) * conv2nm
		# rms_input_all = np.std((self.inst.phase_ncpa + \
		# 						self.inst.phase_wv)[self.inst.aperture==1]) * conv2nm
		rms_res_qs = np.std(res_ncpa_qs[self.inst.aperture>=0.5]) * conv2nm
		rms_res_all = np.std(res_ncpa_all[self.inst.aperture>=0.5]) * conv2nm

		if self.cfg.params.psi_correction_mode is not 'all':
			tmp, _ = self._projectOnModalBasis(res_ncpa_qs)
			rms_res_qs_filt = np.std(tmp[self.inst.aperture>=0.5]) * conv2nm
			tmp, _ = self._projectOnModalBasis(res_ncpa_all)
			rms_res_all_filt = np.std(tmp[self.inst.aperture>=0.5]) * conv2nm
		else:
			rms_res_qs_filt = rms_res_qs
			rms_res_all_filt = rms_res_all

		tmp, _ = self._projectOnModalBasis(self.inst.phase_wv)
		rms_wv = np.std(tmp[self.inst.aperture>=0.5]) * conv2nm

		tmp, _ = self._projectOnModalBasis(self.inst.phase_ncpa_correction)
		rms_corr = np.std(tmp[self.inst.aperture>=0.5]) * conv2nm

		tmp, _ = self._projectOnModalBasis(res_static_ncpa_qs)
		rms_res_static_NCPA_filt = np.std(tmp[self.inst.aperture>=0.5]) * conv2nm

		if verbose:
			self.logger.info('#{0} : Res [QS, QS+WV] = [{1:.0f}, {2:.0f}]'.\
				format(self.iter, rms_res_qs, rms_res_all))
			self.logger.info('#{0} : Res. filt. [QS, QS+WV] = [{1:.0f}, {2:.0f}]'.\
				format(self.iter, rms_res_qs_filt, rms_res_all_filt))
			self.logger.info('#{0} : input WV_f rms  = {1:.0f}'.format(self.iter, rms_wv))
			self.logger.info('#{0} : PSI correction rms = {1:.0f}'.format(self.iter, rms_corr))
			self.logger.info('#{0} : Long-term (static) residual rms = {1:.0f}'.format(self.iter, rms_res_static_NCPA_filt))


		loop_stat = [self.iter]
		loop_stat.append(rms_res_all_filt)
		loop_stat.append(rms_res_qs_filt)
		loop_stat.append(rms_res_all)
		loop_stat.append(rms_res_qs)
		self._loop_stats.append(loop_stat)

#
# if __name__ == '__main__':
# 	# config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_metis_compass.py'
# 	config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_demo_metis_compass.py'
#
# 	psi_sensor = PsiSensor(config_file)
#
# 	psi_sensor.setup()
# 	# Test: doing one iteration
# 	psi_sensor.logger.info('Inputs:')
# 	psi_sensor.evaluateSensorEstimate()
# 	psi_sensor.ncpa_scaling = 1e-3
# 	# psi_sensor.next()
# 	# psi_sensor.evaluateSensorEstimate()
# 	# psi_sensor.next()
# 	# psi_sensor.evaluateSensorEstimate()
# 	# for i in range(10):
# 	# 	psi_sensor.next()
# 	# 	psi_sensor.evaluateSensorEstimate()
# 	psi_sensor.loop()
