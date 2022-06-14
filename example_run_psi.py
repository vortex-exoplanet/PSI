from psi.psiSensor import PsiSensor

# config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_metis_compass.py'
config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_demo_metis_compass.py'

psi_sensor = PsiSensor(config_file)

psi_sensor.setup()
# Test: doing one iteration
psi_sensor.logger.info('Inputs:')
psi_sensor.evaluateSensorEstimate()
# psi_sensor.ncpa_scaling = 1e-3
# psi_sensor.next()
# psi_sensor.evaluateSensorEstimate()
# psi_sensor.next()
# psi_sensor.evaluateSensorEstimate()
# for i in range(10):
# 	psi_sensor.next()
# 	psi_sensor.evaluateSensorEstimate()
psi_sensor.loop()
