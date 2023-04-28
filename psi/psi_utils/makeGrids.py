import numpy as np
from scipy.signal import convolve2d
from .psi_utils import gauss_2Dalt
from hcipy.aperture import make_circular_aperture, make_obstructed_circular_aperture
from hcipy.coronagraphy import VortexCoronagraph
from hcipy import inverse_tikhonov
from hcipy.mode_basis import ModeBasis, make_zernike_basis
from hcipy.optics import Apodizer, make_gaussian_influence_functions, OpticalSystem


def makeFilters(grid, type="back_prop", sigma=0.05, cent_obs=0.27, outer_vin=1.0,
                ravel=True, lD=15):
    """Packages required:
    Circular and make_obstructed_circular_aperture
    gauss_2Dalt
    convolve2d

    HISTORY:
     - 24/03/2022: adding a 'ravel' keywork to return the 2D array. + some re-writing for clarification

    COMMENTS:
     - cent_obs is not used if type=='back_pro' -- default parameter should be None. And warning should be issued. Same applies to outer_vin
     - parameters should be exposed for all 3 types of filters

    TODO: improvement in makeFilters:
    - clean up the makeFilters
    - expose the size of the aperture (the 15lbda/D by default)
    - option between Gaussian filter and simple circular filter

    """
    if type == "back_prop":
        filter_ = make_circular_aperture(lD)(grid)  # 15 is an arbitrary number set by Emiel.
    elif type == "ncpa":
        filter_ = make_obstructed_circular_aperture(1*0.7, 1.8*cent_obs)(grid)
    elif type == "reset":
        filter_ = make_obstructed_circular_aperture(1*outer_vin, 1.8*cent_obs)(grid)
    else:
        print("Type of filter not given!!")
        quit()
    filter_2D = filter_.shaped
    sigma_ = int(sigma*filter_2D.shape[0])
    # sigma_ = int(sigma*np.sqrt(filter_.shape))
    if sigma_ != 0:
        # kernel_ = gauss_2Dalt(size_x=int(np.sqrt(filter_.shape)), sigma_x=sigma_)
        # filter_.shape = (int(np.sqrt(filter_.shape)),int(np.sqrt(filter_.shape)))  # -> why? what ?!
        # filter_ = convolve2d(filter_, kernel_, mode='same')
        kernel_ = gauss_2Dalt(size_x=filter_2D.shape[0], sigma_x=sigma_)
        filter_2D = convolve2d(filter_2D, kernel_, mode='same')
        if ravel:
            filter_1D = filter_2D.ravel()
            return filter_1D
        else:
            return filter_2D
    # return filter_


def makeMatrices(grid_, ao_acts_, aperture_, rcond):
    '''
    HISTORY:
    2022-03-24 (GOX):
            - removing the orthoganlization of ao_modes_ -> this generate strange influence functions
            - renaming 'reoncstruction_normalisation_' to 'rcond'

    TODO: change name to 'makeWfsCalibrationMatrices'
    '''
    ao_modes_ = make_gaussian_influence_functions(
        grid_, ao_acts_, 1.0 / ao_acts_, crosstalk=0.15)  # 1.2						# Create an object containing all the available DM pistons, 1.0 to
    # ao_modes_ = ModeBasis([mode * aperture_ for mode in ao_modes_]).orthogonalized
    ao_modes_ = ModeBasis([mode * aperture_ for mode in ao_modes_])
    transformation_matrix_ = ao_modes_.transformation_matrix
    reconstruction_matrix_ = inverse_tikhonov(transformation_matrix_, rcond)
    return transformation_matrix_, reconstruction_matrix_


def makeOpticalSystem(grid_):
    '''
    TODO: Improve makeOpticalSystem
            - expose Lyot stop parameters
            - option: use the VectorVortex
            - option: RAVC
            - option: APP
    '''
    lyot_stop_ = make_obstructed_circular_aperture(0.98, 0.3)(grid_)
    coro_ = OpticalSystem([VortexCoronagraph(grid_, 2), Apodizer(lyot_stop_)])
    return coro_


def makeZerns(n_zerns, grid, start):
    '''
    TODO: Improve makeZerns and change name to makeModalBasis (?)
            - use orthogonalization on specific aperture + proper normalization
            - different basis: disk harmonics, Karhunen-Loeve, ... ?
            - expose the regularization parameter
            - possibility to use custom basis (given via cube fits) ?
    '''
    zernike_modes_ = make_zernike_basis(n_zerns, 1, grid, start).orthogonalized
    reconstructor_zernike_ = inverse_tikhonov(zernike_modes_.transformation_matrix, 1e-3)
    modal_means_ = np.ones(n_zerns)
    return zernike_modes_, reconstructor_zernike_, modal_means_
