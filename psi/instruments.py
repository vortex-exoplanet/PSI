import os
import abc
from types import SimpleNamespace
import sys
# sys.path.append('')
import hcipy
# sys.path.append('/Users/orban/Projects/METIS/4.PSI/psi_github/')
import psi.psi_utils as psi_utils
from psi.psi_utils.psi_utils import crop_img, resize_img
import numpy as np
from .helperFunctions import LazyLogger
import astropy.io.fits as fits


class GenericInstrument():
    '''
    Abstract Generic instrument

    Define all the necessary methods and properties that are used by
    ``PsiSensor``.
    '''
    def __init__(self, conf):
        '''
        Compute the pupilGrid and the focalGrid within the HCIPy framework.
        Compute phase buffer as HCIPy.Field.
        Define the propagator from pupil to focal plane as Fraunhofer propagator
        '''
        self._size_pupil_grid = conf.npupil
        self._focal_grid_resolution = conf.det_res
        self._focal_grid_size = conf.det_size

        self.pupilGrid = hcipy.make_pupil_grid(self._size_pupil_grid)
        self.focalGrid = hcipy.make_focal_grid(self._focal_grid_resolution,
                                               self._focal_grid_size)

        self._prop = hcipy.FraunhoferPropagator(self.pupilGrid, self.focalGrid)


        self.phase_ncpa = hcipy.Field(0.0, self.pupilGrid)         # knowledge only in Simulation
        self.phase_wv = hcipy.Field(0.0, self.pupilGrid)           # knowledge only in Simulation
        self.phase_ncpa_correction = hcipy.Field(0.0, self.pupilGrid)  # NCPA correction applied

        pass

    @property
    def optical_model(self):
        '''
        HCIPy.OpticalSystem object
            a linear path of optical elements that propagates the wavefront
            forward and backward.
        '''
        return self._optical_model

    @optical_model.setter
    def optical_model(self, model):
        '''
        HCIPy.OpticalSystem object
            a linear path of optical elements that propagate the wavefront
            forward and backward.
        '''
        self._optical_model = model

    @property
    def aperture(self):
        '''
        Entrance pupil of the instrument
        '''
        return self._aperture

    @aperture.setter
    def aperture(self, aper):
        '''
        Setting the entrance pupil of the instrument
        '''
        self._aperture = aper

    @abc.abstractmethod
    def grabWfsTelemetry(self):
        '''
        Grab wavefront sensor telemetry and returns wavefront buffer
        '''
        pass

    @abc.abstractmethod
    def grabScienceImages(self):
        '''
        Grab science images and returns science images buffer
        '''
        pass

    @abc.abstractmethod
    def setNcpaCorrection(self):
        '''
            Set the NPCA phase correction.
        '''
        pass

    @abc.abstractmethod
    def synchronizeBuffers(self):
        '''
        Synchronize science and wfs telemetry buffers
        '''
        pass

    @abc.abstractmethod
    def getNumberOfPhotons(self):
        '''
        Provides an estimate of the total number of photons at
        the entrance pupil plane

        (calibration task)
        '''
        pass


class CompassSimInstrument(GenericInstrument):
    def __init__(self, conf, logger=LazyLogger('CompassInstrument')):
        super().__init__(conf)

        self.logger = logger

        if type(conf) == dict:
            conf = SimpleNamespace(**conf)

        self._setup(conf)

        start_idx = 2011
        self._zero_time_ms = start_idx
        self._current_time_ms = start_idx  # starting time with COMPASS phase_screens
        self._start_time_wfs = start_idx
        self._end_time_wfs = start_idx
        self._start_time_sci_buffer = start_idx
        self._end_time_sci_buffer = start_idx
        self._start_time_last_sci_dit = start_idx
        self._end_time_last_sci_dit = start_idx

    def _setup(self, conf):
        self.wfs_exptime = 1 / conf.ao_framerate
        self.ao_frame_decimation = conf.ao_frame_decimation
        self.sci_exptime = conf.dit

        self.nb_ao_per_sci = conf.nb_ao_frames_per_science
        # self.ncpa_map = conf.ncpa_map
        # self.add_water_vapour = conf.add_water_vapour

        self.wavelength = conf.wavelength
        self._prefix_rp = conf.turb_prefix_rp
        self._prefix_wf = conf.turb_prefix_wf
        self._suffix = conf.turb_suffix
        self._input_folder = conf.turb_folder

        # Aperture definition -- GOX: see also modification later in _setup to cope with slight mismatch with the NCPA maps
        self.aperture = psi_utils.make_COMPASS_aperture(conf.f_aperture,
                                                  npupil=self._size_pupil_grid,
                                                  rot90=True,
                                                  binary=True)(self.pupilGrid)
        # self.aperture = np.rot90(self.aperture)
        self._inst_mode = conf.inst_mode  # which type of imaging system
        if self._inst_mode == 'CVC' or self._inst_mode == 'RAVC':
            self._vc_charge = conf.vc_charge
            self._vc_vector = conf.vc_vector
            # Lyot stop mask definition ...
            self.lyot_stop_mask = psi_utils.make_COMPASS_aperture(conf.f_lyot_stop,
                                                            npupil=self._size_pupil_grid,
                                                            rot90=True)(self.pupilGrid)
            # self.lyot_stop_mask = np.rot90(self.lyot_stop_mask)
        if self._inst_mode == 'RAVC' or self._inst_mode == 'APP':
            self.pupil_apodizer = psi_utils.make_COMPASS_aperture(conf.f_apodizer,
                                                            npupil=self._size_pupil_grid,
                                                            rot90=True)(self.pupilGrid)

        self.noise = conf.noise
        # if self.noise == 1:
        #     pass
        if self.noise == 2:
            self.bckg_level = conf.num_photons_bkg
        self.num_photons = conf.num_photons

        # by default include residual turbulence phase screens
        self.include_residual_turbulence = True
        self.phase_residual = 0

        self.ncpa_dynamic = conf.ncpa_dynamic
        if self.ncpa_dynamic:
            self.ncpa_sampling = conf.ncpa_sampling
        self._input_folder_ncpa = conf.ncpa_folder
        self._prefix_ncpa = conf.ncpa_prefix
        self.ncpa_scaling = conf.ncpa_scaling
        self._initialize_dynamic_ncpa()
        # --- Customizing the entrance aperture due to some small mismatch in the NCPA definition ---
        # nb: see the 'hardcoded' +6 pixels -- same as in psi_utils.loadNCPA
        self.logger.warn('Customizing entrance aperture as a function of NCPA map definition')
        ncpa_file = self._prefix_ncpa + str(self._ncpa_index) + '.fits'
        mask_pupil = fits.getdata(self._input_folder_ncpa + ncpa_file)
        mask_pupil[mask_pupil != 0 ]=1
        # mask_pupil = resize_img(mask_pupil, self._size_pupil_grid)
        # mask_pupil= np.rot90(mask_pupil)
        size_ = self._size_pupil_grid
        mask_pupil = psi_utils.psi_utils.process_screen(mask_pupil, size_+6,
                                                  self.aperture, rotate=True, ncpa_=True)
        mask_pupil = psi_utils.psi_utils.crop_img(mask_pupil, (size_, size_))


        mask_pupil = np.ravel(mask_pupil)
        self.aperture = self.aperture * mask_pupil
        #-------------#

        self.include_water_vapour = conf.wv
        if self.include_water_vapour:
            self.wv_folder = conf.wv_folder
            self.wv_cubename = conf.wv_cubename
            self.wv_sampling = conf.wv_sampling
            self.wv_scaling = conf.wv_scaling
            self._initialize_water_vapour()
        else:
            pass
            # self.phase_wv = 0  # already initialilze in GenericInstrument
        # self.aperture = conf.aperture

        # self.phase_ncpa_correction = 0
        # ....

        # self.start_time = 2011
        # COMPASS units are µm; HCIPy needs rad
        self.conv2rad = 1e3 * (2 * np.pi / self.wavelength * 1e-9)

        # HEEPS cube are in meters:

        self.toto_scaling = 1

    def build_optical_model(self):
        '''
            Building an optical model in HCIPy depending on the instrument mode selected.

            Instrument modes are
                - CVC :  Classical Vector Coronagraph
                - ELT : normal imaging (no coronagraph, no Lyot stop)
                - RAVC : Ring-Apodized Vector Coronagraph
                - APP : Apodized Phase Plate coronagraph

        '''
        if self._inst_mode == 'CVC':
            self.logger.info('Building a Classical Vortex Coronagraph optical model in HCIPy')

            assert self._vc_charge == 2 or self._vc_charge == 4

            if self._vc_vector:
                self._vvc_element = hcipy.VectorVortexCoronagraph(self._vc_charge)
            else:
                self._vvc_element = hcipy.VortexCoronagraph(self.pupilGrid, self._vc_charge)

            self._lyot_stop_element = hcipy.Apodizer(self.lyot_stop_mask)

            self.optical_model = hcipy.OpticalSystem([self._vvc_element,
                                                      self._lyot_stop_element,
                                                      self._prop])
        elif self._inst_mode == 'ELT':
            self.logger.info('Building a simple imager in HCIPy')
            self.optical_model = hcipy.OpticalSystem([self._prop])

        elif self._inst_mode == 'RAVC':

            self.logger.info('Building a Ring-Apodizer Vortex Coronagraph optical model in HCIPy')

            self._ring_apodizer = hcipy.Apodizer(self.pupil_apodizer)

            assert self._vc_charge == 2 or self._vc_charge == 4

            if self._vc_vector:
                self._vvc_element = hcipy.VectorVortexCoronagraph(self._vc_charge)
            else:
                self._vvc_element = hcipy.VortexCoronagraph(self.pupilGrid, self._vc_charge)

            self._lyot_stop_element = hcipy.Apodizer(self.lyot_stop_mask)

            self.optical_model = hcipy.OpticalSystem([self._ring_apodizer,
                                                      self._vvc_element,
                                                      self._lyot_stop_element,
                                                      self._prop])

        elif self._inst_mode == 'APP':
            self.logger.warning('APP not supported')

        # # lyot_stop_mask = hcipy.make_obstructed_circular_aperture(0.98, 0.3)(pupil_grid)
        # # lyot_stop_mask = hp.evaluate_supersampled(hp.circular_aperture(0.95), pupil_grid, 4)
        # lyot_stop_mask = hp.circular_aperture(0.95)
        # lyot_stop = hcipy.Apodizer(lyot_stop_mask)



    def _initialize_dynamic_ncpa(self):
        # -- NCPA should be part of the instrument ... not here ----
        self._ncpa_index = 0
        ncpa_file = self._prefix_ncpa + str(self._ncpa_index) + '.fits'
        size_pupil_grid = int(self.pupilGrid.shape[0])
        self.phase_ncpa = psi_utils.loadNCPA(self.aperture, size_pupil_grid,
                                       file_=ncpa_file,
                                       folder_=self._input_folder_ncpa)
        self.phase_ncpa *= self.ncpa_scaling
        # # compute min max for plot
        # ncpa_min = - np.ptp(self.phase_ncpa) / 2
        # ncpa_max = np.ptp(self.phase_ncpa) / 2

    def _update_dynamic_ncpa(self, current_time):
        '''read/compute a new NCPA map'''

        if (((current_time - self._zero_time_ms)/1e3) % self.ncpa_sampling) == 0:
            self.logger.info('Updating NCPA map')
            ncpa_file = self._prefix_ncpa+str(self._ncpa_index) + '.fits'
            size_pupil_grid = int(self.pupilGrid.shape[0])
            self.phase_ncpa = psi_utils.loadNCPA(self.aperture,
                                           size_pupil_grid,
                                           file_=ncpa_file,
                                           folder_=self._input_folder_ncpa)
            self.phase_ncpa *= self.ncpa_scaling
            self._ncpa_index += 1

    def _initialize_water_vapour(self):
        self._wv_index = 0
        # HEEPS cubes are in meters
        self.conv2rad_wv = (2 * np.pi / self.wavelength)
        self.phase_wv_cube = fits.getdata(self.wv_folder + self.wv_cubename)
        size_pupil_grid = int(self.pupilGrid.shape[0])
        self.phase_wv  = self.conv2rad_wv * \
            psi_utils.process_screen(self.phase_wv_cube[0],
                               size_pupil_grid,
                               self.aperture, rotate=True)
        # folder_wv = '/Users/orban/Projects/METIS/4.PSI/legacy_TestArea/WaterVapour/phases/'
        # file_wv = "cube_Cbasic_20210504_600s_100ms_0piston_meters_scao_only_285_WVLonly_qacits.fits"
        # wave_vapour_cube = fits.getdata(os.path.join(folder_wv, file_wv)) * \
            # 2 * np.pi / wavelength  #* 1e3 * 1e-6
        # pass

    def _update_water_vapour(self, current_time):
        '''read/compute a new NCPA map'''

        if ((current_time - self._zero_time_ms) %  self.wv_sampling) == 0:
            # self.logger.info('Updating WV map')
            size_pupil_grid = int(self.pupilGrid.shape[0])
            self.phase_wv  = self.conv2rad_wv * \
                psi_utils.process_screen(self.phase_wv_cube[self._wv_index],
                                   size_pupil_grid,
                                   self.aperture, rotate=True)
            self._wv_index += 1



    def grabScienceImages(self, nbOfPastSeconds):
        '''
            Grab a buffer of science images
        '''
        self.nbOfSciImages = int(nbOfPastSeconds / self.sci_exptime)
        assert self.nbOfSciImages <= nbOfPastSeconds / self.sci_exptime
        if not(np.isclose(self.nbOfSciImages, nbOfPastSeconds/self.sci_exptime)):
            self.logger.warn('Requested buffer duration is not an integer number of Science DIT')

        nx, ny = self.focalGrid.shape
        image_buffer = np.zeros((self.nbOfSciImages, nx, ny))

        self._start_time_sci_buffer = np.copy(self._current_time_ms)
        # re-initialize timer of single dit
        self._start_time_last_sci_dit = np.copy(self._start_time_sci_buffer)
        self._end_time_last_sci_dit = np.copy(self._start_time_sci_buffer)
        for i in range(self.nbOfSciImages):
            image_buffer[i] = self._grabOneScienceImage()

        self._end_time_sci_buffer = np.copy(self._end_time_last_sci_dit)
        return image_buffer


    def _grabOneScienceImage(self):
        '''
            Compute a single science image: consist of several realisation of the
            residual turbulence (+ NPCA, WV, NCPA_correction)
        '''
        # conversion_COMPASSToNm = 1e3
        # conv = conversion_COMPASSToNm * (2 * np.pi / self.wavelength * 1e-9)
        # conv = 2 * np.pi
        nbOfFrames = int(self.sci_exptime / (self.wfs_exptime * self.ao_frame_decimation))
        deltaTime = (self.wfs_exptime * self.ao_frame_decimation) * 1e3
        timeIdxInMs = np.arange(nbOfFrames) * deltaTime

        # self._start_time_sci = np.copy(self._current_time_ms)
        # file_indices = [str(int(self._current_time_ms + timeIdxInMs[i]))
        #                 for i in range(len(timeIdxInMs))]
        self._start_time_last_sci_dit= np.copy(self._end_time_last_sci_dit)
        file_indices = [str(int(self._start_time_last_sci_dit + timeIdxInMs[i]))
                        for i in range(len(timeIdxInMs))]

        # phase in radians
        file_wf = self._prefix_rp + '_' + file_indices[0] + self._suffix
        phase_pupil = fits.getdata(os.path.join(self._input_folder, file_wf)) * self.conv2rad

        # Remove piston
        phase_pupil = psi_utils.remove_piston(phase_pupil, self.aperture.shaped)
        # conversion to HCIPy
        residual_phase = hcipy.Field(phase_pupil.ravel(), self.pupilGrid)
        wf_post_ = hcipy.Wavefront(np.exp(1j * residual_phase) * self.aperture,
                                   self.wavelength)
        # Setting number of photons
        # wf_post_.total_power = self.num_photons
        # Propagation through the instrument
        efield_fp = self.optical_model(wf_post_)
        img_one = efield_fp.power

        nx, ny = img_one.shaped.shape
        image_cube = np.zeros((nbOfFrames, nx, ny))
        ss = residual_phase.shape[0]
        total_phase_cube = np.zeros((nbOfFrames, ss))

        for i in range(len(file_indices)):
            # self._current_time_ms = self._current_time_ms + timeIdxInMs[i]

            file_wf = self._prefix_rp + '_' + file_indices[i] + self._suffix

            #
            if self.include_residual_turbulence:
                self.phase_residual = fits.getdata(os.path.join(self._input_folder,
                                                           file_wf)) * self.conv2rad

                self.phase_residual = psi_utils.remove_piston(self.phase_residual, self.aperture.shaped)
                self.phase_residual *= self.toto_scaling

            # Update water vapour phase
            if self.include_water_vapour :
                self._update_water_vapour(self._start_time_last_sci_dit + timeIdxInMs[i])

            # Update NCPA phase
            if self.ncpa_dynamic :
                self._update_dynamic_ncpa(self._start_time_last_sci_dit + timeIdxInMs[i])

            # Get current NCPA correction

            total_phase_cube[i] = self.phase_residual.ravel() + \
                self.phase_wv + self.phase_ncpa + self.phase_ncpa_correction

        # Forward propagation and calculation of the image for a sequence of phases
        total_phase_cube = hcipy.Field(total_phase_cube, self.pupilGrid)
        # wf_post_ = hcipy.Wavefront(np.exp(1j * total_phase_cube) * self.aperture, 1)
        wf_post_ = hcipy.Wavefront(np.exp(1j * total_phase_cube) * self.aperture)

        # Setting number of photons
        # ToDo
        wf_post_.total_power = self.num_photons * nbOfFrames
        # Propagation through the instrument
        # TODO: expose 'prop' and 'coro'

        self._image_cube = self.optical_model(wf_post_).power.shaped
        # if vvc:
        # 	image_cube[i] = prop(coro(wf_post_)).power.shaped
        # else:
        # 	image_cube[i] = prop((wf_post_)).power.shaped
        assert len(self._image_cube.shape) == 3
        image = self._image_cube.mean(0)
        # Photometry -- TBC
        if self.noise == 0:
            noisy_image = image
        elif self.noise == 1:
            noisy_image = hcipy.large_poisson(image)
        elif self.noise == 2:
            background_noise = hcipy.large_poisson(self.bckg_level + image*0) - \
                self.bckg_level
            noisy_image = hcipy.large_poisson(image) + background_noise
        # +	np.random.poisson(nb_photons, image.shape)

        self._end_time_last_sci_dit = self._start_time_last_sci_dit + timeIdxInMs[-1] + deltaTime

        return noisy_image


    def grabWfsTelemetry(self, nbOfPastSeconds):
        '''
        Returns
                phase cube in units of radian
        '''
        # self._compass_start_time=2011 # COMPASS 0 indexing in msec

        # conversion_COMPASSToNm = 1e3
        # conv = conversion_COMPASSToNm * (2 * np.pi / self.wavelength * 1e-9)
        nbOfFrames = int(nbOfPastSeconds / (self.wfs_exptime * self.ao_frame_decimation))
        deltaTime = (self.wfs_exptime * self.ao_frame_decimation) * 1e3
        timeIdxInMs = np.arange(nbOfFrames) * deltaTime


        self._start_time_wfs = np.copy(self._current_time_ms)
        file_indices = [str(int(self._current_time_ms + timeIdxInMs[i]))
                        for i in range(len(timeIdxInMs))]

        fname = self._prefix_wf + '_' + file_indices[0] + self._suffix
        phase_pupil = fits.getdata(os.path.join(self._input_folder, fname)) *\
            self.conv2rad

        phase_cube = np.zeros((nbOfFrames, phase_pupil.shape[0], phase_pupil.shape[1]))

        for i in range(len(file_indices)):
            # self._current_time_ms = self._current_time_ms + timeIdxInMs[i]
            file_wf = self._prefix_wf + '_' + file_indices[i] + self._suffix

            # read file
            phase = fits.getdata(os.path.join(self._input_folder, file_wf)) *\
                self.conv2rad
            # remove piston
            phase = psi_utils.remove_piston(phase, self.aperture.shaped)
            phase *= self.toto_scaling

            phase_cube[i] = np.copy(phase)

        self._end_time_wfs = self._current_time_ms + timeIdxInMs[-1] + deltaTime
        return phase_cube

    def setNcpaCorrection(self, phase):
        self.phase_ncpa_correction = self.phase_ncpa_correction + phase

    def synchronizeBuffers(self, wfs_telemetry_buffer, sci_image_buffer):
        '''
            Note:
                wfs_telemetry_buffer & sci_image_buffer are actually not used here.
                To be realistic, one could correlate the tip-tilt in both to sync them.
        '''
        if self._start_time_wfs != self._start_time_sci_buffer:
            self.logger.warn('Start buffers not sync')
            self.logger.debug('Start WFS buffer is {0}'.format(self._start_time_wfs))
            self.logger.debug('Start SCI buffer is {0}'.format(self._start_time_sci_buffer))
            # return 0
        if self._end_time_wfs != self._end_time_sci_buffer:
            self.logger.warn('End buffers not sync')
            self.logger.debug('End WFS buffer is {0}'.format(self._end_time_wfs))
            self.logger.debug('End SCI buffer is {0}'.format(self._end_time_sci_buffer))
            # return 0

        self._current_time_ms = np.copy(self._end_time_wfs)

        # For each science image, calculate a start and stop index for
        #   the wfs telemetry buffer
        telemetry_indexing = [(i * self.nb_ao_per_sci, (i+1) * self.nb_ao_per_sci)
                           for i in range(self.nbOfSciImages)]

        return telemetry_indexing

    def getNumberOfPhotons(self):
        return self.num_photons

class DemoCompassSimInstrument(CompassSimInstrument):

    def __init__(self, conf, logger=LazyLogger('CompassInstrument')):
        super().__init__(conf)

    def _grabOneScienceImage(self):
        '''
            Compute a single science image: consist of several realisation of the
            residual turbulence (+ NPCA, WV, NCPA_correction)
        '''
        # conversion_COMPASSToNm = 1e3
        # conv = conversion_COMPASSToNm * (2 * np.pi / self.wavelength * 1e-9)
        # conv = 2 * np.pi
        nbOfFrames = int(self.sci_exptime / (self.wfs_exptime * self.ao_frame_decimation))
        deltaTime = (self.wfs_exptime * self.ao_frame_decimation) * 1e3
        timeIdxInMs = np.arange(nbOfFrames) * deltaTime

        # self._start_time_sci = np.copy(self._current_time_ms)
        # file_indices = [str(int(self._current_time_ms + timeIdxInMs[i]))
        #                 for i in range(len(timeIdxInMs))]
        self._start_time_last_sci_dit= np.copy(self._end_time_last_sci_dit)
        file_indices = [str(int(self._start_time_last_sci_dit + timeIdxInMs[i]))
                        for i in range(len(timeIdxInMs))]

        # phase in radians
        file_wf = self._prefix_rp + '_'  + self._suffix
        phase_pupil = fits.getdata(os.path.join(self._input_folder, file_wf)) * self.conv2rad

        # Remove piston
        phase_pupil = psi_utils.remove_piston(phase_pupil, self.aperture.shaped)
        # conversion to HCIPy
        residual_phase = hcipy.Field(phase_pupil.ravel(), self.pupilGrid)
        wf_post_ = hcipy.Wavefront(np.exp(1j * residual_phase) * self.aperture,
                                   self.wavelength)
        # Setting number of photons
        # wf_post_.total_power = self.num_photons
        # Propagation through the instrument
        efield_fp = self.optical_model(wf_post_)
        img_one = efield_fp.power

        nx, ny = img_one.shaped.shape
        image_cube = np.zeros((nbOfFrames, nx, ny))
        ss = residual_phase.shape[0]
        total_phase_cube = np.zeros((nbOfFrames, ss))

        for i in range(len(file_indices)):
            # self._current_time_ms = self._current_time_ms + timeIdxInMs[i]

            # file_wf = self._prefix_rp + '_' + file_indices[i] + self._suffix

            #
            if self.include_residual_turbulence:
                self.phase_residual = phase_pupil


            # Get current NCPA correction

            total_phase_cube[i] = self.phase_residual.ravel() + \
                self.phase_wv + self.phase_ncpa + self.phase_ncpa_correction

        # Forward propagation and calculation of the image for a sequence of phases
        total_phase_cube = hcipy.Field(total_phase_cube, self.pupilGrid)
        # wf_post_ = hcipy.Wavefront(np.exp(1j * total_phase_cube) * self.aperture, 1)
        wf_post_ = hcipy.Wavefront(np.exp(1j * total_phase_cube) * self.aperture)

        # Setting number of photons
        # ToDo
        wf_post_.total_power = self.num_photons * nbOfFrames
        # Propagation through the instrument
        # TODO: expose 'prop' and 'coro'

        self._image_cube = self.optical_model(wf_post_).power.shaped
        # if vvc:
        # 	image_cube[i] = prop(coro(wf_post_)).power.shaped
        # else:
        # 	image_cube[i] = prop((wf_post_)).power.shaped
        assert len(self._image_cube.shape) == 3
        image = self._image_cube.mean(0)
        # Photometry -- TBC
        if self.noise == 0:
            noisy_image = image
        elif self.noise == 1:
            noisy_image = hcipy.large_poisson(image)
        elif self.noise == 2:
            background_noise = hcipy.large_poisson(self.bckg_level + image*0) - \
                self.bckg_level
            noisy_image = hcipy.large_poisson(image) + background_noise
        # +	np.random.poisson(nb_photons, image.shape)

        self._end_time_last_sci_dit = self._start_time_last_sci_dit + timeIdxInMs[-1] + deltaTime

        return noisy_image

    def grabWfsTelemetry(self, nbOfPastSeconds):
        '''
        Returns
                phase cube in units of radian
        '''
        # self._compass_start_time=2011 # COMPASS 0 indexing in msec

        # conversion_COMPASSToNm = 1e3
        # conv = conversion_COMPASSToNm * (2 * np.pi / self.wavelength * 1e-9)
        nbOfFrames = int(nbOfPastSeconds / (self.wfs_exptime * self.ao_frame_decimation))
        deltaTime = (self.wfs_exptime * self.ao_frame_decimation) * 1e3
        timeIdxInMs = np.arange(nbOfFrames) * deltaTime


        self._start_time_wfs = np.copy(self._current_time_ms)
        file_indices = [str(int(self._current_time_ms + timeIdxInMs[i]))
                        for i in range(len(timeIdxInMs))]

        fname = self._prefix_wf + '_'  + self._suffix
        phase_pupil = fits.getdata(os.path.join(self._input_folder, fname)) *\
            self.conv2rad

        phase_cube = np.zeros((nbOfFrames, phase_pupil.shape[0], phase_pupil.shape[1]))

        for i in range(len(file_indices)):
            # self._current_time_ms = self._current_time_ms + timeIdxInMs[i]
            # file_wf = self._prefix_wf + '_' + file_indices[i] + self._suffix

            # # read file
            # phase = fits.getdata(os.path.join(self._input_folder, file_wf)) *\
            #     self.conv2rad
            # # remove piston
            # phase = psi_utils.remove_piston(phase, self.aperture.shaped)
            # phase *= self.toto_scaling

            phase_cube[i] = np.copy(phase_pupil)

        self._end_time_wfs = self._current_time_ms + timeIdxInMs[-1] + deltaTime
        return phase_cube

class HcipySimInstrument(GenericInstrument):
    pass


class ErisInterfaceOffline(GenericInstrument):
    pass


# if __name__ == '__main__':
#     from configParser import loadConfiguration
#     config_file = '/Users/orban/Projects/METIS/4.PSI/psi_github/config/config_metis_compass.py'
#     cfg = loadConfiguration(config_file)
#     inst = CompassSimInstrument(cfg.params)
#     inst.build_optical_model()
