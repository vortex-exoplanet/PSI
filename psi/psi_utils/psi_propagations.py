import numpy as np
from hcipy.optics import Wavefront
from hcipy.propagation import Propagator
from hcipy import large_poisson


def prop_image(wf_post_coro_, ncpa_,
               correction_, prop_, coro_,
               wv=0, noise_=True, n_pho=8.88e+04/1000, background_subtract=False):
    '''
    HISTORY:
    - [GOX] wf_post_coro = wf_post_coro_.copy() and then use wf_post_coro. Otherwise the 'external' field get modified...

    TODO: check photometry (/formula) which is not clear to me

    '''

    #
    wf_post_coro = wf_post_coro_.copy()
    # * np.exp(1j * wv) # Add NCPA and correction
    wf_post_coro.electric_field *= np.exp(1j * ncpa_) * np.exp(1j * correction_)
    if coro_ is None:
        img_oneEx_ = prop_(wf_post_coro).power
    else:
        img_oneEx_ = prop_(coro_(wf_post_coro)).power
    if noise_ == True:
        image_ = large_poisson(img_oneEx_) + np.random.poisson(n_pho, img_oneEx_.shape)
    else:
        image_ = img_oneEx
    if background_subtract == True:
        image_ -= np.full(image_.shape, n_pho)
        image_ = np.abs(image_)
    # image_.shape = (int(np.sqrt(image_.shape)),int(np.sqrt(image_.shape)))
    # image_ = np.roll(image_,1,1)
    # image_ = image_.ravel()
    return image_


def prop_wf(wf_post_ao_, aperture_, prop_, coro_, recon_m, trans_m):
    '''
    HISTORY
    - [GOX] remove subtraction by aperture_
                    (create artificial phase wrapping and not apparent reason for it)
    - [GOX] adding again the aperture_ subtraction -> required for good behavior...
    '''
    #
    # * aperture is different between METIS and ERIS
    wfs_measurement_ = recon_m.dot(
        np.angle(wf_post_ao_.electric_field / wf_post_ao_.electric_field.mean()) * aperture_)
    reconstructed_pupil = aperture_ * np.exp(1j * trans_m.dot(wfs_measurement_))
    reconstructed_pupil /= np.exp(1j * np.angle(reconstructed_pupil.mean()))
    reconstructed_pupil -= aperture_

    if coro_ is None:
        reconstructed_electric_field_ = prop_(Wavefront(reconstructed_pupil)).electric_field
    else:
        reconstructed_electric_field_ = prop_(coro_(Wavefront(reconstructed_pupil))).electric_field
    return reconstructed_electric_field_
