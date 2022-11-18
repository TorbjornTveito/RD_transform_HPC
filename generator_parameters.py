R_M = 1737.4e3

#radar parameters
c = 299_792_458


class radar_paramaters():
    def __init__(self, obs_frequency, radar_lon, radar_lat, radar_el, baud_length, gain, tx_power, noise_temp, duty_cycle, obsdur, obsdate, make_SAR, do_area, do_SPRCS, do_SLAW, RD_decimation_factor = 1):
        self.RD_decimation_factor = RD_decimation_factor                # decimation factor for generating low-resolution datasets. should be integer. 1 if not used
        self.obs_frequency = obs_frequency                              # central frequency of transmit pulse IN HERTZ
        self.wavelength = c / obs_frequency                             # central wavelength of transmit pulse IN METERS
        self.radar_lon = radar_lon                                      # longitude of radar
        self.radar_lat = radar_lat                                      # latitude of radar
        self.radar_el = radar_el                                        # elevation of radar in meters
        self.baud_length = baud_length                                  # baud length of transmitted pulse in seconds
        self.range_res = RD_decimation_factor * c * baud_length / 2     # range resolution in meters
        self.gain = gain                                                # gain of radar - does not currently support different transmitters and receivers
        self.tx_power = tx_power                                        # transmit power of radar transmitter
        self.noise_temp = noise_temp                                    # noise temperature of receiver in kelvins
        self.duty_cycle = duty_cycle                                    # duty cycle of radar transmitter as a percentage
        self.obsdur = obsdur                                          # duration of observation in seconds
        self.obsdate = obsdate                                          # Time of observation
        self.make_SAR = make_SAR                                        # whether to do SAR calculations
        self.do_area = do_area                                          # whether to include the effects of SAR pixel area
        self.do_SPRCS = do_SPRCS                                        # include signal-per-radar-crossection in SAR calculations
        self.do_SLAW = do_SLAW                                          # calculate scattering law effects


list_of_axes =[
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 11, 1, 0), True, False, False, False, 1),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 13, 0, 0), True, False, False, False, 1),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 18, 2, 0), True, False, False, False, 1),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 19, 5, 0), True, False, False, False, 1),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 20, 7, 0), True, False, False, False, 1),
]

decimated_list_of_axes = [
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 11, 1, 0), True, False, False, False, 512),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 13, 0, 0), True, False, False, False, 512),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 18, 2, 0), True, False, False, False, 512),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 19, 5, 0), True, False, False, False, 512),
    radar_paramaters(930e6, -71.4882, 42.6233, 130.0, 1/1e6, 3600, 2e6, 80, 0.125, 3600, (2022, 3, 20, 7, 0), True, False, False, False, 512),
]




#observation parameters

#job parameters