def planck_function ( self , wavelength , temperature ) :
    """ Calculate spectral radiance from Planck ’s law .
    Args :
    wavelength ( float or array ): Wavelength in meters
    temperature ( float ): Temperature in Kelvin
    
    Returns :
    float or array : Spectral radiance in W.m -2.sr -1.um -1
    """



def planck_inv ( self , wavelength , spectral_radiance ) :
    """ Calculate temperature from spectral radiance using Planck ’s law .
    
    Args :
    wavelength ( float ): Wavelength in meters
    spectral_radiance ( float ): Spectral radiance in W.m -2.sr -1.um -1
    
    Returns :
    float : Temperature in Kelvin
    """

    
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










