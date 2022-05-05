import numpy as np
import astropy.io.fits as fits
import os
# from ..field import CartesianGrid, UnstructuredCoords, make_hexagonal_grid, Field
# from .generic import *

# These two lines maybe needed in the future for making custom apertures so I have left them in
from hcipy.field import CartesianGrid, UnstructuredCoords, make_hexagonal_grid, Field
from hcipy.aperture.generic import *
from .psi_utils import crop_img, resize_img

def make_vlt_aperture():
    pass


def make_subaru_aperture():
    pass


def make_lbt_aperture():
    pass


def make_elt_aperture():
    pass


def make_COMPASS_aperture(fname,
                          npupil=256,
                          # input_folder='/Users/matt/Documents/METIS/TestArea/fepsi/COMPASSPhaseScreens/Test/',
                          # file_name='mask_256.fits',
                          nimg=720):
    '''Create an aperture from a COMPASS product.

    Parameters
    ----------
    npupil : scalar
            Number of pixels across each dimension of the array.
    input_folder : string
            Location of the aperture file from COMPASS.
    nimg : scalar
            The size the aperture file needs to be cut down to as it comes with some padding. 720 should be the default unless something
            changes in the COMPASS products.

    Returns
    -------
    Field generator
            The resized COMPASS aperture.
    '''

    # mask = fits.getdata(os.path.join(input_folder, 'mask_256.fits'))
    mask = fits.getdata(fname)
    if mask.shape[0] < nimg:
        mask = crop_img(mask, nimg, verbose=False)
    mask_pupil = resize_img(mask, npupil)
    #mask_pupil[mask_pupil<0.8] = 0
    # mask_pupil = mask_pupil.transpose() # Testing for wind direction dependencies. Should be commented out.
    aperture = np.ravel(mask_pupil)

    # nimg = 720
    # npupil = 256

    def func(grid):
        return Field(aperture, grid)
    return func


def make_ERIS_aperture(npupil=256,
                       input_folder='',
                       file_name_='VLT_ERIS_pupil_D_tightly_undersized_ZEMAX_corrected_SkyBaffle.fits',
                       nimg=2500):
    '''Create an aperture from a COMPASS product.

    Parameters
    ----------
    npupil : scalar
            Number of pixels across each dimension of the array.
    input_folder : string
            Location of the aperture file from COMPASS.
    nimg : scalar
            The size the aperture file needs to be cut down to as it comes with some padding. 720 should be the default unless something
            changes in the COMPASS products.

    Returns
    -------
    Field generator
            The resized COMPASS aperture.
    '''

    mask = fits.getdata(os.path.join(input_folder, file_name_))
    if mask.shape[0] < nimg:
        mask = crop_img(mask, nimg, verbose=False)
    mask_pupil = resize_img(mask, npupil)
    #mask_pupil[mask_pupil<0.8] = 0
    # mask_pupil = mask_pupil.transpose() # Testing for wind direction dependencies. Should be commented out.
    aperture = np.ravel(mask_pupil)

    # nimg = 720
    # npupil = 256

    def func(grid):
        return Field(aperture, grid)
    return func
