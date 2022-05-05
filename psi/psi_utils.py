
from astropy.io import fits
import numpy as np
import os
from skimage.transform import resize
import warnings
from hcipy.optics import Wavefront
from hcipy.field import Field


def crop_img(img, new_size, margin=0, verbose=False):
    ''' Crop an img to a new size. Handles even and odd sizes.
    Can add an optional margin of length 1, 2 (x,y) or 4 (x1,x2,y1,y2).
    '''
    requirement = "new_size must be an int or a tuple/list of size 2."
    assert type(new_size) in [int, tuple, list], requirement
    if type(new_size) is int:
        (x1, y1) = (new_size, new_size)
    else:
        assert len(new_size) == 2, requirement
        (x1, y1) = new_size
    (x2, y2) = img.shape
    if not np.any(np.array([x1, y1]) < np.array([x2, y2])):
        if verbose == True:
            print('crop size is larger than img size')
    else:
        # determine cropping region
        dx = int((x2 - x1)/2)
        dy = int((y2 - y1)/2)
        cropx = (dx, dx) if (x2-x1) % 2 == 0 else (dx+1, dx)
        cropy = (dy, dy) if (y2-y1) % 2 == 0 else (dy+1, dy)
        # check for margins
        requirement2 = "margin must be an int or a tuple/list of size 2 or 4."
        assert type(margin) in [int, tuple, list], requirement2
        if type(margin) is int:
            (mx1, mx2, my1, my2) = (margin, margin, margin, margin)
        elif len(margin) == 2:
            (mx1, mx2, my1, my2) = (margin[0], margin[0], margin[1], margin[1])
        else:
            assert len(margin) == 4, requirement2
            (mx1, mx2, my1, my2) = margin
        # crop image
        img = img[cropx[0]-mx1:-cropx[1]+mx2, cropy[0]-my1:-cropy[1]+my2]
    return img


def resize_img(img, new_size, preserve_range=True, mode='reflect',
               anti_aliasing=True):
    ''' Resize an image. Handles even and odd sizes.
    '''
    requirement = "new_size must be an int or a tuple/list of size 2."
    assert type(new_size) in [int, tuple, list], requirement
    if type(new_size) is int:
        new_size = (new_size, new_size)
    else:
        assert len(new_size) == 2, requirement
    assert img.ndim in [2, 3], 'image must be a frame (2D) or a cube (3D)'
    if img.ndim == 3:
        new_size = (len(img), *new_size)
    if new_size != img.shape:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # when anti_aliasing=False, and NANs
            img = np.float32(resize(np.float32(img), new_size,
                                    preserve_range=preserve_range, mode=mode, anti_aliasing=anti_aliasing))
    return img


def pad_img(img, padded_size, pad_value=0):
    ''' Pad an img with a value (default is zero). Handles even and odd sizes.
    '''
    requirement = "padded_size must be an int or a tuple/list of size 2."
    assert type(padded_size) in [int, tuple, list], requirement
    if type(padded_size) is int:
        (x1, y1) = (padded_size, padded_size)
    else:
        assert len(padded_size) == 2, requirement
        (x1, y1) = padded_size
    (x2, y2) = img.shape
    # determine padding region
    assert not (x1 < x2 or y1 < y2), "padding region can't be smaller than image size."
    dx = int((x1 - x2)/2)
    dy = int((y1 - y2)/2)
    padx = (dx, dx) if (x1-x2) % 2 == 0 else (dx+1, dx)
    pady = (dy, dy) if (y1-y2) % 2 == 0 else (dy+1, dy)
    # pad image
    img = np.pad(img, [padx, pady], mode='constant', constant_values=pad_value)
    return img


def process_screen(phase_, grid_size_, aperture_, rotate=False, ncpa_=False):

    if rotate == True and ncpa_ == False:
        return resize_img(np.rot90(phase_), grid_size_).ravel() * aperture_
    elif rotate == True and ncpa_ == True:
        return resize_img(np.rot90(phase_), grid_size_)
    elif rotate == False and ncpa_ == True:
        return resize_img(phase_, grid_size_)
    else:
        return resize_img(phase_, grid_size_).ravel() * aperture_


def load_file(folder, file, grid_size, aperture_):
    '''
    TODO: specific to COMPASS/METIS ?
    '''
    ncpa_ = fits.getdata(os.path.join(folder, file)) * 2 * np.pi / \
        wavelength * 1e3 * 1e-6  # mm * rad/um * m/mm * um/m
    # The NCPA from Tibor come in 512 so need to be resized.
    ncpa_ = process_screen(ncpa_, grid_size_, aperture_, rotate=True)
    return ncpa_


def loadNCPA(aperture_,
             size_,
             file_="DIFF_rep_1_field_0.fits",
             folder_="/Users/matt/Documents/METIS/TestArea/fepsi/NCPA_Tibor/",
             wavelength_=3e-6):
    """The NCPA from Tibor come in 512 so need to be resized. Their aperture is slightly smaller than
    COMPASS, 8pxs to the edge after resizing rather than 6px. +6 is then for enlarging a bit further and
    then use crop to bring it back down to the right size.

    GOX comments on units:
    -----------------------
    HCIPy needs radians
    Tibor phase screens are expressed in [mm]
    -> hence load_NCPA includes a (1e3 1e-6) = 1e-3 to convert millimeters to meters
    """
    ncpa_ = fits.getdata(os.path.join(folder_, file_)) * 2 * np.pi / \
        wavelength_ * 1e3 * 1e-6  # mm * rad/um * m/mm * um/m
    ncpa_ = process_screen(ncpa_, size_+6, aperture_, rotate=True, ncpa_=True)
    ncpa_ = crop_img(ncpa_, (size_, size_))
    ncpa_ = ncpa_.ravel() * aperture_
    return ncpa_


def load_and_process(file, in_fold, mask_pup, grid, pow=1e3, wavelength=1.0, piston=True):
    '''
    PARAMETERS
    -----------
    file: str
        fits file name of the phase screen to load.
        units should be in µm (TBC)
    in_fold : str
        path name where the phase screen is
    mask_pup : Field
        telescope  aperture defined as a HCIPy Field
    grid: CartesianGrid
        HCIPy pupil grid
    pow : float
        setting the total power (amplitude normalisation) --> why is it needed???
    wavelength: float
        wavelength for the conversion from µm to rad (TBC)
    piston: bool
        remove the piston in the phase screen over the aperture. Default is True

    GOX modifications:
    -------------------
        - Adding this doc
        - 'if remove_piston' -> 'if piston'
        - fixing mask_pup -> mask_pup.shaped as input for the piston removal
        - multiplying the phase arg by the mask_pup

    - TODO: check wavelength normalisation
    - TODO: check utility of setting the total_power
    '''
    ### Load phase residuals ###
    phase_pupil = fits.getdata(os.path.join(in_fold, file)) * 2 * np.pi / wavelength
    # phase_pupil  = phase_pupil.transpose() # Testing for wind direction dependencies. Should be commented out.
    ### Transform the phase residuals into electric field and then into the right format ##
    if piston == True:
        phase_pupil = remove_piston(phase_pupil, mask_pup.shaped)
    phase_residual_ao = Field(phase_pupil.ravel(), grid)
#     wf_post_ = Wavefront(np.exp(1j * phase_residual_ao*mask_pup) * mask_pup, wavelength)
    wf_post_ = Wavefront(np.exp(1j * phase_residual_ao*mask_pup) * mask_pup)

    wf_post_.total_power = pow
    return wf_post_


def get_contrast_curve(frame, radius_resolution):
    arr_dim = int(np.sqrt(frame.shape[0]))
    x = np.linspace(0, arr_dim, arr_dim) - arr_dim/2
    y = np.linspace(0, arr_dim, arr_dim) - arr_dim/2
    xv, yv = np.meshgrid(x, y)
    rd = np.sqrt(xv**2 + yv**2)
    #radii = np.arange(0,np.max(rd)+radius_resolution, radius_resolution)
    radii = np.arange(0, arr_dim, radius_resolution)
    mean_flux = np.empty(0)
    for ring_ind in np.arange(np.size(radii)-1):
        rad_inner = radii[ring_ind]
        rad_outer = radii[ring_ind+1]
        rad_in = rd > rad_inner
        rad_out = rd < rad_outer
        rad = rad_in*rad_out
        mean_flux = np.append(mean_flux, np.mean(frame[np.ravel(rad)]))
    return mean_flux


def remove_piston(data, mask_):
    data += -np.mean(data[mask_ != 0.0])  # Remove piston
    data[mask_ == 0.0] = 0.0
    return data


def gauss_2Dalt(size_x=400, size_y=None, sigma_x=5, sigma_y=None, x0=None, y0=None, theta=0):
    # 2D Guassian which can be elongated in x,y and angled.
    if size_y == None:
        size_y = size_x
    if sigma_y == None:
        sigma_y = sigma_x
    if x0 == None:
        x0 = size_x//2
    if y0 == None:
        y0 = x0
    size_x = int(size_x)
    size_y = int(size_y)
    x_g, y_g = np.meshgrid(np.arange(0, size_x), np.arange(0, size_y))
    x_g = x_g - x0
    y_g = y_g - y0
    a = np.cos(theta)**2 / (2*sigma_x**2) + np.sin(theta)**2 / (2*sigma_y**2)
    b = -np.sin(2*theta) / (4*sigma_x**2) + np.sin(2*theta) / (4*sigma_y**2)
    c = np.sin(theta)**2 / (2*sigma_x**2) + np.cos(theta)**2 / (2*sigma_y**2)
    exp_part = a*x_g**2 - 2*b*x_g*y_g + c*y_g**2
    return 1/(2*np.pi*sigma_x*sigma_y) * np.exp(-exp_part)


def resetVariables(a):
    return np.copy(a), 0, 0, 0+0j


def resetPSIVariables():
    return 0, 0, 0, 0
