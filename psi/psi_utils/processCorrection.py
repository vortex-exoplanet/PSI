import numpy as np

def processCorrection(ncpa_estimate_, correction_, correction_normalisation,
					  aperture_, ncpa_filter_aperture_, recon_zern_recon, \
					  zern_recon, modal_means_init_, modal_means_,
					  reset=False, resetFilter=None, ZernCorrection=True, ms=0):

	m_corrections      = recon_zern_recon.dot(ncpa_estimate_ * aperture_)

	# Lets have a look at the zernikes which are produced. Starting from the 4th zernike (no piston or tip/tilt)
	if np.size(ms) != 1:
		m_mean = np.mean(m_corrections)
		m_std  = np.std(m_corrections)
		ms += m_corrections#*modal_means_

	ncpa_estimate_temp = np.copy(ncpa_estimate_)
	new_correction_max = np.max(zern_recon.transformation_matrix.dot(m_corrections*modal_means_init_))

	lo_freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
	lo_freq_corr      = np.max(lo_freq_estimate) / new_correction_max
	lo_freq_estimate *= lo_freq_corr
	hi_freq_estimate  = (ncpa_estimate_ - lo_freq_estimate) * ncpa_filter_aperture_

	# Correcting by power rather than max
	# ncpa_estimate_temp = np.copy(ncpa_estimate_)
	# new_correction_sum = np.sum(zern_recon.transformation_matrix.dot(m_corrections*modal_means_init_))
	# lo_freq_corr      = np.sum(lo_freq_estimate) / new_correction_sum
	# lo_freq_estimate *= lo_freq_corr

	if ZernCorrection==True:
		new_mask = aperture_-ncpa_filter_aperture_	# Define a mask of just the problem areas. Let everything else be high frequency
		new_estimate = (lo_freq_estimate*new_mask*0.5 + ncpa_estimate_*ncpa_filter_aperture_)
		new_estimate *= np.max(new_estimate) / new_correction_max
		correction_ += -correction_normalisation * new_estimate

		# correction_ += -(correction_normalisation * new_mask * lo_freq_estimate * 0.5) -(correction_normalisation * hi_freq_estimate)
	else:
		# freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
		# freq_estimate *= np.max(freq_estimate) / new_correction_max
		# correction_   += -(correction_normalisation * (ncpa_estimate_-freq_estimate))
		# freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
		# freq_estimate *= lo_freq_corr
		# hi_freq_estimate = (ncpa_estimate_ - freq_estimate) * (aperture_-ncpa_filter_aperture_)
		# # freq_estimate *= np.max(freq_estimate) / new_correction_max
		# correction_   += -(correction_normalisation * freq_estimate) +(correction_normalisation * hi_freq_estimate)/10.

		m_corrections      = recon_zern_recon.dot(ncpa_estimate_ * aperture_)
		ncpa_estimate_temp = np.copy(ncpa_estimate_)
		new_correction_max = np.max(zern_recon.transformation_matrix.dot(m_corrections*modal_means_init_))

		lo_freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
		lo_freq_estimate *= np.max(lo_freq_estimate) / new_correction_max
		hi_freq_estimate  = (ncpa_estimate_ - lo_freq_estimate)/10. * ncpa_filter_aperture_

		correction_ += -(correction_normalisation * lo_freq_estimate)# -(correction_normalisation * hi_freq_estimate)

		# correction_   += -(correction_normalisation * lo_freq_estimate) -(correction_normalisation * hi_freq_estimate)

	if reset == True:
		try:
			correction_ *= resetFilter
		except ValueError:
			print("No reset filter given.")

	# correction_ += -correction_normalisation*ncpa_estimate_
	if np.size(ms) != 1:
		return correction_, ms
	else:
		return correction_

def processCorrection_Original(ncpa_estimate_, correction_,
							   correction_normalisation, aperture_,
							   ncpa_filter_aperture_, recon_zern_recon, 
							   zern_recon, modal_means_init_, modal_means_,
							   reset=False, resetFilter=None,
							   ZernCorrection=True, ms=0):

	m_corrections      = recon_zern_recon.dot(ncpa_estimate_ * aperture_)
	ncpa_estimate_temp = np.copy(ncpa_estimate_)
	new_correction_max = np.max(zern_recon.transformation_matrix.dot(m_corrections*modal_means_init_))

	if np.size(ms) != 1:
		m_mean = np.mean(m_corrections)
		m_std  = np.std(m_corrections)
		ms += m_corrections#*modal_means_

	lo_freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
	lo_freq_estimate *= np.max(lo_freq_estimate) / new_correction_max
	hi_freq_estimate  = (ncpa_estimate_ - lo_freq_estimate)/10. * ncpa_filter_aperture_

	correction_ += -(correction_normalisation * lo_freq_estimate) -(correction_normalisation * hi_freq_estimate)

	# correction_ += -correction_normalisation*ncpa_estimate_
	if np.size(ms) != 1:
		return correction_, ms
	else:
		return correction_

#####
# Original

# m_corrections      = recon_zern_recon.dot(ncpa_estimate_ * aperture_)
# ncpa_estimate_temp = np.copy(ncpa_estimate_)
# new_correction_max = np.max(zern_recon.transformation_matrix.dot(m_corrections*modal_means_init_))
#
# lo_freq_estimate  = zern_recon.transformation_matrix.dot(m_corrections*modal_means_)
# lo_freq_estimate *= np.max(lo_freq_estimate) / new_correction_max
# hi_freq_estimate  = (ncpa_estimate_ - lo_freq_estimate)/10. * ncpa_filter_aperture_
#
# correction_ += -(correction_normalisation * lo_freq_estimate) -(correction_normalisation * hi_freq_estimate)
