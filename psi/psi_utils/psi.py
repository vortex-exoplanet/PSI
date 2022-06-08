import numpy as np
from hcipy.optics import Wavefront


def PSI(ref_, img_,
        I_sum, phi_sum, phi_2, phi_I,
        filter_, grid, prop_, i_,
        returnEpup=False):
    '''
            HISTORY:
             - 16/03/2022 : adding the possibility to return the full electric field estimate
             - 13/04/2022 : addin case i==0
    '''
    #
    phi_I += ref_ * img_
    phi_2 += np.abs(ref_)**2
    phi_sum += ref_
    I_sum += img_

    # TODO: psi algebra and backward propagation is in principle
    #       only needed once ensemble average is calculated
    #       TBC (if it holds true, that would relax the computation)
    if i_ == 0:
        psi_estimate_ = phi_I / phi_2 ** 2
    else:
        psi_estimate_ = (phi_I - phi_sum * I_sum / (i_+1)) / (phi_2 - np.abs(phi_sum / (i_+1))**2)

    wf = Wavefront(psi_estimate_ * filter_)

    # This line solves the phase wrapping
    wf.electric_field *= np.exp(-2j * grid.as_('polar').theta)
    # TODO: in coronagraphic mode,
    #   this should be the backward propagation over the full coronagaraph
    #   and not sur a backwar Fraunhofer prop.
    pup = prop_.backward(wf)
    ncpa_estimate_ = pup.electric_field.imag

    if returnEpup:
        return ncpa_estimate_, I_sum, phi_sum, phi_2, phi_I, psi_estimate_, pup
    else:
        return ncpa_estimate_, I_sum, phi_sum, phi_2, phi_I, psi_estimate_
