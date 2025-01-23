import numpy as np
from scipy.integrate import quad

class Radiation:
    def __init__(self):
        # Constants
        self.h = 6.62607015e-34  # Planck constant in J.s
        self.c = 3.0e8  # Speed of light in m/s
        self.k = 1.380649e-23  # Boltzmann constant in J/K


    
    def planck_function ( self , wavelength , temperature ) :
        """ Calculate spectral radiance from Planck ’s law .
        Args :
        wavelength ( float or array ): Wavelength in meters
        temperature ( float ): Temperature in Kelvin
        
        Returns :
        float or array : Spectral radiance in W.m -2.sr -1.um -1
        """
        
        # Planck function
        spectral_radiance = (2.0 * self.h * self.c**2) / (wavelength**5) * (1.0 / (np.exp((self.h * self.c) 
                                                                                          / (wavelength * self.k * temperature)) - 1.0))
        spectral_radiance_um = spectral_radiance * 1e-6
        
        return spectral_radiance_um
    
    
    
    def planck_inv ( self , wavelength , spectral_radiance ) :
        """ Calculate temperature from spectral radiance using Planck ’s law .
        
        Args :
        wavelength ( float ): Wavelength in meters
        spectral_radiance ( float ): Spectral radiance in W.m -2.sr -1.um -1
        
        Returns :
        float : Temperature in Kelvin
        """
        # Convert spectral radiance to W.m^-2.sr^-1.m^-1
        spectral_radiance_m = spectral_radiance * 1e6

        # Use Planck function inverse formula to calculate temperature
        temperature = (self.h * self.c) / (wavelength * self.k * np.log((2.0 * self.h * self.c**2) / (wavelength**5 * spectral_radiance_m) + 1.0))

        return temperature
    
        
    def brightness_temperature ( self , radiance , band_center , band_width ) :
        """ Calculate brightness temperature for a given rectangular bandpass .
        
        Implements numerical integration over a rectangular bandpass defined
        by its center wavelength and width . The bandpass is assumed to have
        unity transmission within its bounds and zero outside .
        
        Args :
        radiance ( float ): Observed radiance in W.m -2.sr -1.um -1
        band_center ( float ): Center wavelength of bandpass in meters
        band_width ( float ): Width of bandpass in meters
        
        Returns :
        float : Brightness temperature in Kelvin
        
        Raises :
        ValueError : If band_width <= 0 or band_center <= band_width /2
        """
        if band_width <= 0 or band_center <= band_width / 2:
            raise ValueError("Invalid band_width or band_center")

        # Define the integration limits
        lower_limit = band_center - band_width / 2
        upper_limit = band_center + band_width / 2

        # Integrate the Planck function over the bandpass
        def integrand(wavelength):
            return self.planck_function(wavelength, radiance)

        integrated_radiance, _ = quad(integrand, lower_limit, upper_limit)

        # Calculate the average radiance over the bandpass
        average_radiance = integrated_radiance / band_width

        # Calculate the brightness temperature
        brightness_temp = self.planck_inv(band_center, average_radiance)

        return brightness_temp
    
    
    
    def radiance ( self , temperature , band_center , band_width ) :
        """ Calculate band - integrated radiance for a given temperature and
        rectangular bandpass .
        
        Integrates Planck function over a rectangular bandpass defined
        by its center wavelength and width . The bandpass is assumed to
        have unity transmission within its bounds and zero outside .
        
        Args :
        temperature ( float ): Temperature in Kelvin
        band_center ( float ): Center wavelength of bandpass in meters
        band_width ( float ): Width of bandpass in meters
        
        Returns :
        float : Band - integrated radiance in W.m -2.sr -1
        
        Raises :
        ValueError : If temperature <= 0 , band_width <= 0 , or
        band_center <= band_width /2
        """
        if temperature <= 0 or band_width <= 0 or band_center <= band_width / 2:
            raise ValueError("Invalid temperature, band_width, or band_center")

        # Define the integration limits
        lower_limit = band_center - band_width / 2
        upper_limit = band_center + band_width / 2

        # Integrate the Planck function over the bandpass
        def integrand(wavelength):
            return self.planck_function(wavelength, temperature)

        integrated_radiance, _ = quad(integrand, lower_limit, upper_limit)

        # Calculate the band-integrated radiance
        band_integrated_radiance = integrated_radiance / band_width

        return band_integrated_radiance
    
    
    def calculate_NEDT ( self , temperature , NER , band_center , band_width ) :
        """ Calculate the noise - equivalent differential temperature ( NEDT )
        for given scene temperature and noise - equivalent radiance (NER ).
        
        Uses numerical derivative of band - integrated radiance with respect
        to temperature to determine the temperature uncertainty corresponding
        to the NER.
        
        Args :
        temperature ( float ): Scene temperature in Kelvin
        NER ( float ): Noise - equivalent radiance in W.m -2. sr -1
        band_center ( float ): Center wavelength of bandpass in meters
        band_width ( float ): Width of bandpass in meters
        
        Returns :
        float : NEDT in Kelvin
        
        Raises :
        ValueError : If temperature <= 0 , NER <= 0 , band_width <= 0 ,
        or band_center <= band_width /2
        """
        if temperature <= 0 or NER <= 0 or band_width <= 0 or band_center <= band_width / 2:
            raise ValueError("Invalid temperature, NER, band_width, or band_center")

        # Calculate the band-integrated radiance at the given temperature
        radiance_at_temp = self.radiance(temperature, band_center, band_width)

        # Define a small temperature increment for numerical differentiation
        delta_temp = 1e-3  # Small temperature change in Kelvin

        # Calculate the band-integrated radiance at temperature + delta_temp
        radiance_at_temp_plus_delta = self.radiance(temperature + delta_temp, band_center, band_width)

        # Calculate the numerical derivative of radiance with respect to temperature
        d_radiance_d_temp = (radiance_at_temp_plus_delta - radiance_at_temp) / delta_temp

        # Calculate the NEDT
        NEDT = NER / d_radiance_d_temp

        return NEDT
    
    
    
    
    




