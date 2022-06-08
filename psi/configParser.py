import traceback
import sys
# sys.path.append('/Users/orban/Projects/METIS/4.PSI/psi_github/')
from .helperFunctions import LazyLogger

class ConfigurationError(Exception):
    pass

class Parameters(object):
    """
    Args:
        filename (string): The name of the configuration file



    Description of the config parameters
        ========================      ===================
        **Required Parameter**        **Description**
        ------------------------      -------------------
        ``npupil``                    int: number of pixels of the pupil
        ``det_size``                  -
        ``det_res``                   -
        ``instrument``                -
        ``inst_mode``                 -
        ``vc_charge``                 -
        ``vc_vector``                 -
        ``f_aperture``                -
        ``f_lyot_stop``               -
        ``dit``                       -
        ``ao_framerate``              -
        ``ao_frame_decimation``       -
        ``psi_framerate``             -
        ``psi_nb_iter``               -
        ``psi_correction_mode``       -
        ``psi_nb_modes``              -
        ``psi_start_mode_idx``        -
        ``ncpa_expected_rms``         -
        ``save_loop_statistics``      -
        ``save_phases_screens``       -
        ``save_basedir``              -
        ========================      ===================

        ======================      ===================
        **Optional Parameter**      **Description**
        ----------------------      -------------------
        ``noise``                   -
        ``mag``                     -
        ``wavelength``              - req ?
        ``flux_zpt``                -
        ``flux_bckg``               -
        ``ncpa_dynamic``            -
        ``ncpa_sampling``           -
        ``ncpa_scaling``            -
        ``ncpa_expected_rms``       -
        ``ncpa_folder``             -
        ``ncpa_prefix``             -
        ``turb_folder``             -
        ``turb_prefix_rp``          -
        ``turb_prefix_wf``          -
        ``turb_suffix``             -
        ``wv``                      -
        ``wv_folder``               -
        ``wv_cube_name``            -
        ``wv_sampling``             -
        ``wv_scaling``              -
        ======================      ===================
    """

    def __init__(self, filename, logger=LazyLogger('Params')):
        self.filename = filename
        self.logger = logger

    def readfile(self):

        #Exec the config file, which should contain a dict ``simConfiguration``
        try:
            with open(self.filename) as file_:
                exec(file_.read(), globals())
        except:
            traceback.print_exc()
            raise ConfigurationError(
                    "Error loading config file: {}".format(self.filename))

        # self.configDict = simConfiguration
        self.params = conf

    def check_parameters(self):
        # check, based on the band, that the zeropoint is according to defined 'constants'
        #  otherwise print a 'warning'

        assert self.params.det_size >= self.params.psi_filt_radius


        # PSI correction mode checks
        if self.params.psi_correction_mode == 'zern':
            if hasattr(self.params, 'psi_nb_modes') is False:
                default_nb_modes = 20
                self.logger.warn('Setting default psi_nb_modes to {0}'.\
                    format(default_nb_modes))
                self.params.psi_nb_modes = default_nb_modes
            if hasattr(self.params, 'psi_start_mode_idx') is False:
                default_start_idx = 4
                self.logger.warn('Setting default psi_start_mode_idx to {0}'.\
                    format(default_start_idx))
                self.params.psi_start_mode_idx = default_start_idx



        # Check if using ``CompassSimInstrument``
        if self.params.instrument == 'CompassSimInstrument':
            if os.path.isfile(self.params.f_aperture) is False:
                self.logger.error('No aperture file, cannot proceed')
                raise ConfigurationError("No aperture file")
        pass

        # Saving
        if self.params.save_phase_screens and not(self.params.save_loop_statistics):
            self.logger.warn('Setting save results to True')
            self.params.save_results = True

    def compute_parameters(self):
        self.params.nb_ao_frames_per_science = int(self.params.dit / \
            (self.params.ao_frame_decimation / self.params.ao_framerate))

        # if Simulation
        self.params.num_photons = self.params.dit * self.params.flux_zpt * \
            10**(-0.4 * self.params.mag)

        self.params.num_photons_bkg = self.params.dit * self.params.flux_bckg



def loadConfiguration(filename):
        cfg_obj = Parameters(filename)

        cfg_obj.readfile()

        cfg_obj.check_parameters()
        cfg_obj.compute_parameters()

        return cfg_obj
