import numpy as np
# from ..field import CartesianGrid, UnstructuredCoords, make_hexagonal_grid, Field
# from .generic import *
from hcipy.field import CartesianGrid, UnstructuredCoords, make_hexagonal_grid, Field	# These two lines maybe needed in the future for making custom apertures so I have left them in
from hcipy.aperture.generic import *

def make_vlt_aperture():
	pass

def make_subaru_aperture():
	pass

def make_lbt_aperture():
	pass

def make_elt_aperture():
	pass

def make_COMPASS_aperture(npupil=256, input_folder='/Users/matt/Documents/METIS/TestArea/fepsi/COMPASSPhaseScreens/Test/', nimg=720):
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

	mask = fits.getdata(os.path.join(input_folder, 'mask_256.fits'))
	if mask.shape[0] < nimg:
		mask = crop_img(mask, nimg, verbose=False)
	mask_pupil = resize_img(mask, npupil)
	#mask_pupil[mask_pupil<0.8] = 0
	#mask_pupil = mask_pupil.transpose() # Testing for wind direction dependencies. Should be commented out.
	aperture = np.ravel(mask_pupil)

	nimg = 720
	npupil = 256

	def func(grid):
		return Field(aperture, grid)
	return func

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
    assert  img.ndim in [2, 3], 'image must be a frame (2D) or a cube (3D)'
    if img.ndim == 3:
        new_size = (len(img), *new_size)
    if new_size != img.shape:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") # when anti_aliasing=False, and NANs
            img = np.float32(resize(np.float32(img), new_size, \
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
    assert not (x1<x2 or y1<y2), "padding region can't be smaller than image size."
    dx = int((x1 - x2)/2)
    dy = int((y1 - y2)/2)
    padx = (dx, dx) if (x1-x2)%2==0 else (dx+1, dx)
    pady = (dy, dy) if (y1-y2)%2==0 else (dy+1, dy)
    # pad image
    img = np.pad(img, [padx, pady], mode='constant', constant_values=pad_value)
    return img

def crop_img(img, new_size, margin=0, verbose=True):
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
    if not np.any(np.array([x1,y1]) < np.array([x2,y2])):
        if verbose == True:
            print('crop size is larger than img size')
    else:
        # determine cropping region
        dx = int((x2 - x1)/2)
        dy = int((y2 - y1)/2)
        cropx = (dx, dx) if (x2-x1)%2==0 else (dx+1, dx)
        cropy = (dy, dy) if (y2-y1)%2==0 else (dy+1, dy)
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
