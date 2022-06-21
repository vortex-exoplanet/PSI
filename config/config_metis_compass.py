


# from heeps.util.download_from_gdrive import extract_zip
import astropy.units as u
import os
import numpy as np
from types import SimpleNamespace
# import proper
# proper.print_it = False


# def read_config(verbose=False, **update_conf):
_tmp_dir = '/Users/orban/Projects/METIS/4.PSI/psi_github/data/'

conf = dict(

    npupil = 256,                        # number of pixels of the pupil
    det_size = 15,                      # [lam/D] size of the detector plane array
    det_res = 4,                         # [px/ (lbda/D)] number of pixels per resolution element

    # --- Which type of instrument to use --
    # Must be a class present in ``instruments.py``
    instrument = 'CompassSimInstrument',

    # =======
    #   Observing Modes
    # =======
    #    0. mode = 'ELT'  for no coronagraph (only telescope)
    #    1. mode = 'CVC'  for Classical Vortex Coronagraph
    #    (2. mode = 'RAVC' for Ring Apodized Vortex Coronagraph)
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
    # # RAVC L-band: ls_RAVC_L_285_dRext=0.0477_dRint=0.02_dRspi=0.0249.fits

    # f_lyot_stop = _tmp_dir + 'pupil/ls_CVC_L_285_dRext=0.0291_dRint=0.08_dRspi=0.0317.fits', # lyot stop file
    f_lyot_stop = _tmp_dir + 'pupil/ls_CVC_L_285_dRext=0.0209_dRint=0.09_dRspi=0.0245.fits',

    # RAVC amptlidue apodization
    f_apodizer = _tmp_dir + 'pupil/apo_ring_r=0.5190_t=0.7909.fits',

    # ======
    #    Photometry
    # ======
    noise = 2  ,                        # 0: no noise, 1: photon noise only, 2: photon noise + background noise
    # add_bckg = False,                   # true means background flux and photon noise are added
    mag = 3,                            # star magnitude at selected band
    # mag_ref = 0,                        # reference magnitude for star and background fluxes
    wavelength = 3.81e-6   ,             # [m] wavelength
    flux_zpt = 8.999e+10,               # [e-/s] zeropoint HCI-L long, mag 0 (Jan 21, 2020)
    flux_bckg = 8.878e+4,              # [e-/s/pix]
    dit = 0.1,                          # [s] science detector integration time

    #bands = 'L', #, 'M', 'N1', 'N2'],

    # [GOX]  this should be somewhere else: this is METIS default value and not supposed to be modified
    #           -> move to a 'constants.py' file or something of the like
    # NB: 'band_specs' is not used by the code. Here for reference
    band_specs = {
        'L': {'lam': 3.81e-6,
            # 'pscale': 5.47,
            'flux_star': 8.999e+10,                 # HCI-L long
            'flux_bckg': 8.878e+04},
        'M': {'lam': 4.79e-6,
            # 'pscale': 5.47,
            'flux_star': 2.452e+10,                 # CO ref
            'flux_bckg': 6.714e+05},
        'N1': {'lam': 8.70e-6,
            # 'pscale': 6.79,
            'flux_star': 3.684e+10,                 # GeoSnap N1
            'flux_bckg': 4.725e+07},
        'N2': {'lam': 11.33e-6,
            # 'pscale': 6.79,
            'flux_star': 3.695e+10,                 # GeoSnap N2
            'flux_bckg': 1.122e+08},
        'N1a': {'lam': 8.67e-6,
            # 'pscale': 10.78,
            'flux_star': 2.979e+10,                 # Aquarius N1
            'flux_bckg': 9.630e+07},
        'N2a': {'lam': 11.21e-6,
            # 'pscale': 10.78,
            'flux_star': 2.823e+10,                 # Aquarius N2
            'flux_bckg': 2.142e+08}
        },

    # ======
    #  AO parameters
    # ======
    ao_framerate = 1000 ,        # [Hz] framerate of the AO loop
    ao_frame_decimation = 10,    # Decimation of the WFS telemetry use by PSI, e.g. if =10, we use 1 every 10 WF frame


    # =========
    #  PSI
    # =========
    psi_framerate = 1,           # [Hz] framerate of the psi correction
    psi_nb_iter = 30,            # number of iterations.

    # How is the PSI estimate process before correction:
    #   1. all     : no projection or filtering
    #   2. zern    : projection on theoretical Zernike mode (circ ap.) and modal control
    #   3. dh    : disk harmonics

    psi_correction_mode = 'zern',
    psi_nb_modes = 100,           # (if modal) nb of modes
    psi_start_mode_idx = 4,        # (if modal) index of first mode. with Zernike, 4 means no piston and tip/tilt

    psi_skip_limit = None,         # [nm rms] value above which the psi_correction will be skipped.
                                  # set to None if no skip limit

    # Focal plane filtering
    psi_filt_sigma = 0.05,
    psi_filt_radius = 10,          # [lbda/D]


    # ============
    #   NCPA
    #       Only in simulation (CompassSimInstrument and HcipySimInstrument)
    # ============
    ncpa_dynamic =  True ,
    ncpa_sampling = 100,             # [s] Dyn NCPA sampling
    ncpa_scaling = 1,               # scaling factor, if want to increase level
    ncpa_expected_rms = 100,        # expected NCPA in [nm]

    ncpa_folder = '/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/NCPA_Tibor/',
    ncpa_prefix = "DIFF_rep_1_field_",  # NB assumes units are in mm
    # ncpa_suffix = '.fits'


    # =============
    #   Residual turbulence
    #       Only in simulation with CompassSimInstrument (offline)
    turb_folder = ('/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/'
                'COMPASSPhaseScreens/ThirdAttempt_Processed/'),
    turb_prefix_rp = 'Residual_phase_screen_',      # NB: assumes units are in µm
    turb_prefix_wf = 'Reconstructed_wavefront_',    # NB: assumes units are in µm
    turb_suffix = 'ms_256.fits',


    # =============
    #   Water vapour seeing
    wv = False,
    wv_folder = '/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/WaterVapour/phases/',
    wv_cubename = 'cube_Cbasic_20210504_600s_100ms_0piston_meters_scao_only_285_WVLonly_qacits.fits',  # NB assume units are in meters
    wv_sampling = 100,      #[ms] sampling of the cube
    wv_scaling = 1,          # scale factor, if want to change the level
    # =============
    # Saving results
    save_loop_statistics = True,
    save_phase_screens = True,
    save_basedir = '/Users/orban/Projects/METIS/4.PSI/psi_results/',

    check_psi_convergence = False,

)
    # sort alphabetically
conf = {k: v for k, v in sorted(conf.items())}
conf = SimpleNamespace(**conf)
# from attrdict import AttrDict
# conf = AttrDict(conf)
    # return conf
