"""Radiative transfer functions for ASTR 5830.

This module provides functions for calculating radiative transfer through
planetary atmospheres, including transmission along arbitrary paths and
cloud scattering calculations.
"""

import numpy as np
#from scipy.interpolate import interp1d
#from radiation_fundamentals import RadiationFundamentals
from fundamentals import Radiation
# Initialize radiation class
rad = Radiation()

# Planck function
planck_function = rad.planck_function

def calc_path_segment_curved(z1: float, 
                           z2: float, 
                           impact_parameter: float,
                           planet_radius: float) -> float:
    """Calculate path length through a spherically symmetric atmospheric layer.
    
    Most basic geometric calculation - path length between two altitudes z1 and z2
    for a given impact parameter (closest approach distance from planet center).
    
    Args:
        z1, z2: Altitudes of segment endpoints in meters
        impact_parameter: Closest approach distance from planet center in meters
        planet_radius: Planet radius in meters
    
    Returns:
        Path length through segment in meters
    """
    r1 = planet_radius + z1
    r2 = planet_radius + z2
    
    # Handle case where path doesn't intersect layer
    if impact_parameter > min(r1, r2):
        return 0.0
        
    # Use law of cosines to find path length
    cos_theta1 = np.sqrt(1 - (impact_parameter/r1)**2)
    cos_theta2 = np.sqrt(1 - (impact_parameter/r2)**2)
    
    return abs(r1*cos_theta1 - r2*cos_theta2)

def integrate_opacity_nadir(opacity_profile: np.ndarray,
                          altitude_grid: np.ndarray,
                          z1: float,
                          z2: float,
                          viewing_angle: float,
                          planet_radius: float) -> float:
    """Calculate integrated opacity along nadir/off-nadir path.
    
    Args:
        opacity_profile: Vertical profile of opacity per unit length (1/m)
        altitude_grid: Altitudes corresponding to opacity profile (m)
        z1: Start altitude of path (m)
        z2: End altitude of path (m)
        viewing_angle: Angle from nadir (radians)
        planet_radius: Planet radius (m)
        
    Returns:
        Integrated optical depth along path
    """
    # Convert viewing angle to local path angle accounting for curvature
    r1 = planet_radius + z1
    r2 = planet_radius + z2
    
    # Use Snell's law in spherical geometry
    sin_alpha = (r1/r2) * np.sin(viewing_angle)
    alpha = np.arcsin(sin_alpha)
    
    # Calculate path length
    dz = z2 - z1
    ds = dz / np.cos(alpha)
    
    # Interpolate opacity to midpoint
    z_mid = (z1 + z2) / 2
    opacity = np.interp(z_mid, altitude_grid, opacity_profile)
    
    # Return optical depth increment
    return opacity * ds

def cloud_thermal_brightness(surface_temperature: float,
                           wavelength: float, 
                           optical_depth: float,
                           single_scatter_albedo: float,
                           emission_angle: float,
                           cloud_temperature: float,
                           cloud_altitude: float,
                           planet_radius: float) -> float:
    """Calculate cloud brightness at thermal wavelengths.
    
    [previous docstring content]
    """
    # Calculate surface and cloud emission using Planck function
    surface_radiance = planck_function(wavelength, surface_temperature)
    cloud_radiance = planck_function(wavelength, cloud_temperature)
    
    # Set up opacity profile for the cloud layer
    dz = 100  # meters
    altitude_grid = np.array([0, cloud_altitude-dz/2, cloud_altitude+dz/2])
    opacity_profile = np.array([0, optical_depth/dz, 0])
    
    # Get slant path optical depth through cloud
    tau = integrate_opacity_nadir(opacity_profile, altitude_grid,
                                cloud_altitude-dz/2, cloud_altitude+dz/2,
                                emission_angle, planet_radius)
    trans = np.exp(-tau)
    
    # Calculate cloud emissivity correctly
    cloud_emissivity = tau * (1 - single_scatter_albedo)
    
    # Components of thermal radiance:
    # 1. Surface emission transmitted through cloud
    surface_component = surface_radiance * trans
    
    # 2. Cloud thermal emission 
    emission_component = cloud_radiance * cloud_emissivity
    
    # 3. Surface emission scattered by cloud
    scatter_component = surface_radiance * tau * single_scatter_albedo
    
    # Total radiance
    return surface_component + emission_component + scatter_component

#------------------------------------------
# Functions to be implemented by students:
#------------------------------------------

def integrate_opacity_limb(opacity_profile: np.ndarray,
                         altitude_grid: np.ndarray,
                         z1: float,
                         z2: float,
                         impact_parameter: float,
                         planet_radius: float) -> float:
    """Calculate integrated opacity along a limb path.
    
    Args:
        opacity_profile: Vertical profile of opacity per unit length (1/m)
        altitude_grid: Altitudes corresponding to opacity profile (m)
        z1: Start altitude of path (m)
        z2: End altitude of path (m) 
        impact_parameter: Closest approach distance from planet center (m)
        planet_radius: Planet radius (m)
    
    Returns:
        Integrated optical depth along limb path
        
    Notes:
        Use calc_path_segment_curved to compute geometric path lengths
    """
     # Calculate path length
    path_length = calc_path_segment_curved(z1, z2, impact_parameter, planet_radius)
    
    # Interpolate opacity to midpoint
    z_mid = (z1 + z2) / 2
    opacity = np.interp(z_mid, altitude_grid, opacity_profile)
    
    # Return optical depth increment
    return opacity * path_length
    
    """

    [alternate (but more complicated) implementation]
    
    # Calculate path segment lengths
    path_segments = np.linspace(z1, z2, num=100)
    path_lengths = np.array([calc_path_segment_curved(z1, z, impact_parameter, planet_radius) for z in path_segments])
    
    # Interpolate opacity profile
    interp_func = interp1d(altitude_grid, opacity_profile, bounds_error=False, fill_value="extrapolate")
    interpolated_opacity = interp_func(path_segments)
    
    # Integrate opacity along the path
    integrated_opacity = np.trapz(interpolated_opacity * path_lengths, path_segments)
    
    return integrated_opacity
    #raise NotImplementedError("Students must implement this function")
    """
    
def surf_transmission(opacity_profile: np.ndarray,
                     angle: float,
                     spacecraft_alt: float,
                     planet_radius: float,
                     altitude_grid: np.ndarray,
                     direction: str = 'from') -> float:
    """Calculate transmission along path to/from surface.
    
    Args:
        opacity_profile: Vertical profile of opacity per unit length (1/m)
        angle: Surface incidence/emission angle (radians)
        spacecraft_alt: Altitude of spacecraft (m)
        planet_radius: Planet radius (m)
        altitude_grid: Altitudes corresponding to opacity profile (m)
        direction: Either 'to' (downward) or 'from' (upward) surface
    
    Returns:
        Total transmission along the path
        
    Raises:
        ValueError: If angle >= pi/2 or direction not 'to'/'from'
        
    Notes:
        Use integrate_opacity_nadir for calculations
    """

    if angle >= np.pi / 2:
        raise ValueError("Angle must be less than π/2 radians.")
    if direction not in ['to', 'from']:
        raise ValueError("Direction must be either 'to' or 'from'.")
    
    if direction == 'to':
        z1 = spacecraft_alt
        z2 = 0.0
    else: #direction == 'from'
        z1 = 0.0
        z2 = spacecraft_alt
    
    # Calculate integrated opacity along the path
    integrated_opacity = integrate_opacity_nadir(opacity_profile, altitude_grid, z1, z2, angle, planet_radius)
    
    # Calculate transmission
    transmission = np.exp(-integrated_opacity)
    
    return transmission
    
    #raise NotImplementedError("Students must implement this function")

def cloud_visible_brightness(surface_albedo: float,
                           solar_flux: float,
                           solar_zenith: float,
                           emission_angle: float, 
                           azimuth_angle: float,
                           optical_depth: float,
                           single_scatter_albedo: float,
                           asymmetry_parameter: float) -> float:
    """Calculate cloud brightness at visible wavelengths using single-scattering.
    
    Implements single-scattering approximation for cloud brightness where
    sunlight is scattered by cloud particles. Includes both direct scattered
    sunlight and surface-reflected light scattered by cloud.
    
    The radiance calculation should implement:
    L = F₀μ₀ϖ₀P(g)/4π × [1 - exp(-τ/μ₀ - τ/μ)] + 
        A_sF₀μ₀/π × exp(-τ/μ₀) × exp(-τ/μ)
    
    Args:
        surface_albedo: Surface Lambert albedo (0-1)
        solar_flux: Incident solar flux at top of atmosphere (W/m2)
        solar_zenith: Solar zenith angle (radians)
        emission_angle: Viewing angle from nadir (radians)
        azimuth_angle: Relative azimuth between sun and viewing direction (radians)
        optical_depth: Cloud optical depth
        single_scatter_albedo: Single scattering albedo (0-1) 
        asymmetry_parameter: Asymmetry parameter g (-1 to 1)
    
    Returns:
        Cloud brightness in W/m2/sr
        
    Notes:
        Use the Henyey-Greenstein phase function:
        P(g) = (1 - g^2)/(1 + g^2 - 2g*cos(theta))^(3/2)
        where theta is the scattering angle
    """
    raise NotImplementedError("Students must implement this function")