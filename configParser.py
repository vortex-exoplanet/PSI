import traceback
from helperFunctions import LazyLogger

class ConfigurationError(Exception):
    pass

class Parameters(object):
    """
    The configuration class holding all simulation configuration information

    This class is used to load the parameter dictionary from file, instantiate each configuration object and calculate some other parameters from the parameters given.

    The configuration file given to this class must contain a python dictionary, named ``simConfiguration``. This must contain other dictionaries for each sub-module of the system, ``Sim``, ``Atmosphere``, ``Telescope``, ``WFS``, ``LGS``, ``DM``, ``Science``. For the final 4 sub-dictionaries, each entry must be formatted as a list (or numpy array) where each value corresponds to that component.

    The number of components on the module will only depend on the number set in the ``Sim`` dict. For example, if ``nGS`` is set to 2 in ``Sim``, then in the ``WFS`` dict, each parameters must have at least 2 entries, e.g. ``subaps : [10,10]``. If the parameter has more than 2 entries, then only the first 2 will be noted and any others discarded.

    Descriptions of the available parameters for each sub-module are given in that that config classes documentation

    Args:
        filename (string): The name of the configuration file

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
