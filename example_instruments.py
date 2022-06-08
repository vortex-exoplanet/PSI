from psi.configParser import loadConfiguration
from psi.instruments import CompassSimInstrument
config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_metis_compass.py'
cfg = loadConfiguration(config_file)
inst = CompassSimInstrument(cfg.params)
inst.build_optical_model()
