


# from heeps.util.download_from_gdrive import extract_zip
import astropy.units as u
import os
import numpy as np
from types import SimpleNamespace
# import proper
# proper.print_it = False


# def read_config(verbose=False, **update_conf):
_tmp_dir = './data/'

conf = dict(

    npupil = 256,                        # number of pixels of the pupil
    det_size = 15,                      # [lam/D] size of the detector plane array
    det_res = 4,                         # [px/ (lbda/D)] number of pixels per resolution element

    # --- Which type of instrument to use --
    # Must be a class present in ``instruments.py``
    instrument = 'DemoCompassSimInstrument',

    # =======
    #   Observing Modes
    # =======
    #    0. mode = 'ELT'  for no coronagraph (only telescope)
    #    1. mode = 'CVC'  for Classical Vortex Coronagraph
    #    2. mode = 'RAVC' for Ring Apodized Vortex Coronagraph
    #    (3. mode = 'APP'  for Apodizing Phase Plate)
    inst_mode = 'ELT',                  # HCI instrument mode
    vc_charge = 2,                      # (CVC and RAVC only) vortex topological charge
    vc_vector = False,                  # (CVC and RAVC only) simulate a vector vortex instead of a scalar one

    # _tmp_dir = '/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/'
    #                  'COMPASSPhaseScreens/Test/'
    #f_aperture = _tmp_dir + 'mask_256.fits',
    f_aperture = _tmp_dir + 'pupil/ELT_fullM1.fits',
    # add parameters to create aperture if no file is given

    # # updates of pupil stops:
    # # CVC L-band: ls_CVC_L_285_dRext=0.0209_dRint=0.09_dRspi=0.0245.fits
    # # CVC N2-band: ls_CVC_N2_119_dRext=0.0268_dRint=0.09_dRspi=0.0357.fits
    # # RAVC L-band: ls_RAVC_L_285_dRext=0.0477_dRint=0.04_dRspi=0.0249.fits

    # RAVC LS L-band:
    f_lyot_stop = _tmp_dir + 'pupil/ls_RAVC_L_285_dRext=0.0477_dRint=0.04_dRspi=0.0249.fits',

    # RAVC amptlidue apodization
    f_apodizer = _tmp_dir + 'pupil/apo_ring_r=0.5190_t=0.7909.fits',


    # ======
    #    Photometry
    # ======
    noise = 2  ,                        # 0: no noise,
                                        # 1: photon noise only,
                                        # 2: photon noise + background noise
    mag = 3,                            # star magnitude at selected band
    wavelength = 3.81e-6   ,            # [m] wavelength
    flux_zpt = 8.999e+10,               # [e-/s] zeropoint HCI-L long, mag 0 (Jan 21, 2020)
    flux_bckg = 8.878e+4,               # [e-/s/pix]
    dit = 0.1,                          # [s] science detector integration time

    # ======
    #  AO parameters
    # ======
    ao_framerate = 1000 ,        # [Hz] framerate of the AO loop
    ao_frame_decimation = 10,    # Decimation of the WFS telemetry use by PSI,
                                 # e.g. if =10, we use 1 every 10 WF frame


    # =========
    #  PSI
    # =========
    psi_framerate = 1,           # [Hz] framerate of the psi correction
    psi_nb_iter = 35,            # number of iterations.

    # How is the PSI estimate process before correction:
    #   1. all     : no projection or filtering
    #   2. zern    : projection on theoretical Zernike mode (circ ap.) and modal control
    #   3. dh    : disk harmonics
    psi_correction_mode = 'zern',
    psi_nb_modes = 100,             # (if modal) nb of modes
    psi_start_mode_idx = 4,         # (if modal) index of first mode. with Zernike,
                                    #  4 means no piston and tip/tilt


    # Focal plane filtering
    psi_filt_sigma = 0.05,
    psi_filt_radius = 15,           # [lbda/D]

    ncpa_expected_rms = 100,        # expected NCPA in [nm]

    # ============
    #   NCPA
    #       Only in simulation (CompassSimInstrument and HcipySimInstrument)
    # ============
    ncpa_dynamic =  False ,
    ncpa_sampling = 100,             # [s] Dyn NCPA sampling
    ncpa_scaling = 1,                # scaling factor, if want to increase level

    ncpa_folder = './data/phase_screens/',
    ncpa_prefix = "DIFF_rep_1_field_",  # NB assumes units are in mm

    # =============
    #   Residual turbulence
    #       Only in simulation with CompassSimInstrument (offline)
    turb_folder = './data/phase_screens/',
    turb_prefix_rp = 'Residual_phase_screen_',      # NB: assumes units are in µm
    turb_prefix_wf = 'Reconstructed_wavefront_',    # NB: assumes units are in µm
    turb_suffix = '2011ms_256.fits',


    # =============
    #   Water vapour seeing
    wv = False,

    # =============
    # Saving results
    save_loop_statistics = True,
    save_phase_screens = True,
    save_basedir = './psi_results/',

)
    # sort alphabetically
conf = {k: v for k, v in sorted(conf.items())}
conf = SimpleNamespace(**conf)
