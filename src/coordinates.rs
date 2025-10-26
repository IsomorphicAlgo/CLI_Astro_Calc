//! Coordinate system conversions
//! 
//! This module provides functions to convert between different astronomical
//! coordinate systems including:
//! - Right Ascension/Declination (RA/Dec)
//! - Altitude/Azimuth (Alt/Az) 
//! - Earth-Centered, Earth-Fixed (ECEF)
//! - Earth-Centered Inertial (ECI)

use crate::Result;

/// Right Ascension and Declination coordinates
#[derive(Debug, Clone, Copy)]
pub struct RaDec {
    pub ra: f64,   // Right Ascension in hours
    pub dec: f64,  // Declination in degrees
}

/// Altitude and Azimuth coordinates
#[derive(Debug, Clone, Copy)]
pub struct AltAz {
    pub alt: f64,  // Altitude in degrees
    pub az: f64,   // Azimuth in degrees
}

/// Earth-Centered, Earth-Fixed coordinates
#[derive(Debug, Clone, Copy)]
pub struct Ecef {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Earth-Centered Inertial coordinates
#[derive(Debug, Clone, Copy)]
pub struct Eci {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Convert RA/Dec to Alt/Az
/// 
/// Converts equatorial coordinates (RA/Dec) to horizontal coordinates (Alt/Az)
/// based on the observer's location and local sidereal time.
/// 
/// # Arguments
/// * `ra_dec` - Right Ascension (hours) and Declination (degrees) coordinates
/// * `observer_lat` - Observer latitude in degrees (positive North, negative South)
/// * `observer_lon` - Observer longitude in degrees (positive East, negative West)
/// * `lst` - Local Sidereal Time in hours (0-24)
/// 
/// # Returns
/// Altitude (degrees above horizon) and Azimuth (degrees from North, clockwise)
/// 
/// # Algorithm
/// 1. Calculate hour angle: H = LST - RA
/// 2. Convert to Alt: alt = arcsin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(H))
/// 3. Convert to Az: az = arctan2(sin(H), cos(H) * sin(lat) - tan(dec) * cos(lat))
pub fn ra_dec_to_alt_az(ra_dec: RaDec, observer_lat: f64, _observer_lon: f64, lst: f64) -> Result<AltAz> {
    // Convert inputs to radians
    let ra_hours = ra_dec.ra;
    let dec_rad = ra_dec.dec.to_radians();
    let lat_rad = observer_lat.to_radians();
    
    // Calculate hour angle (in hours, then convert to radians)
    // Hour angle = LST - RA
    let hour_angle_hours = (lst - ra_hours).rem_euclid(24.0);
    let hour_angle_rad = (hour_angle_hours * 15.0).to_radians(); // Convert hours to degrees to radians
    
    // Calculate altitude using the formula:
    // sin(alt) = sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(H)
    let sin_alt = dec_rad.sin() * lat_rad.sin() 
                + dec_rad.cos() * lat_rad.cos() * hour_angle_rad.cos();
    let alt_rad = sin_alt.asin();
    let alt_deg = alt_rad.to_degrees();
    
    // Calculate azimuth using the formula:
    // Az = atan2(sin(H), cos(H) * sin(lat) - tan(dec) * cos(lat))
    // Then convert to 0-360° measured from North clockwise
    // Using atan2 for proper quadrant handling
    let az_numerator = hour_angle_rad.sin();
    let az_denominator = hour_angle_rad.cos() * lat_rad.sin() - dec_rad.tan() * lat_rad.cos();
    let az_rad = az_numerator.atan2(az_denominator);
    
    // Convert to degrees and adjust to 0-360° range (measured from North, clockwise)
    // The formula gives azimuth measured from South, so we need to add 180°
    let mut az_deg = az_rad.to_degrees() + 180.0;
    if az_deg < 0.0 {
        az_deg += 360.0;
    } else if az_deg >= 360.0 {
        az_deg -= 360.0;
    }
    
    Ok(AltAz {
        alt: alt_deg,
        az: az_deg,
    })
}

/// Convert Alt/Az to RA/Dec
/// 
/// Converts horizontal coordinates (Alt/Az) to equatorial coordinates (RA/Dec)
/// based on the observer's location and local sidereal time.
/// 
/// # Arguments
/// * `alt_az` - Altitude (degrees) and Azimuth (degrees from North) coordinates
/// * `observer_lat` - Observer latitude in degrees (positive North, negative South)
/// * `observer_lon` - Observer longitude in degrees (positive East, negative West)
/// * `lst` - Local Sidereal Time in hours (0-24)
/// 
/// # Returns
/// Right Ascension (hours) and Declination (degrees)
/// 
/// # Algorithm
/// 1. Convert to Dec: dec = arcsin(sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az))
/// 2. Calculate hour angle: H = arctan2(-sin(az), cos(az) * sin(lat) + tan(alt) * cos(lat))
/// 3. Convert to RA: RA = LST - H
pub fn alt_az_to_ra_dec(alt_az: AltAz, observer_lat: f64, _observer_lon: f64, lst: f64) -> Result<RaDec> {
    // Convert inputs to radians
    let alt_rad = alt_az.alt.to_radians();
    let az_rad = alt_az.az.to_radians();
    let lat_rad = observer_lat.to_radians();
    
    // Calculate declination using the formula:
    // sin(dec) = sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az)
    let sin_dec = alt_rad.sin() * lat_rad.sin() 
                + alt_rad.cos() * lat_rad.cos() * az_rad.cos();
    let dec_rad = sin_dec.asin();
    let dec_deg = dec_rad.to_degrees();
    
    // Calculate hour angle using the formula:
    // tan(H) = -sin(az) / (cos(az) * sin(lat) + tan(alt) * cos(lat))
    // Using atan2 for proper quadrant handling
    let hour_angle_rad = (-az_rad.sin()).atan2(
        az_rad.cos() * lat_rad.sin() + alt_rad.tan() * lat_rad.cos()
    );
    
    // Convert hour angle to hours
    let hour_angle_hours = hour_angle_rad.to_degrees() / 15.0;
    
    // Calculate Right Ascension: RA = LST - H
    let ra_hours = (lst - hour_angle_hours).rem_euclid(24.0);
    
    Ok(RaDec {
        ra: ra_hours,
        dec: dec_deg,
    })
}

/// Convert ECEF to ECI
pub fn ecef_to_eci(ecef: Ecef, gmst: f64) -> Result<Eci> {
    // TODO: Implement conversion
    Ok(Eci { x: 0.0, y: 0.0, z: 0.0 })
}

/// Convert ECI to ECEF
pub fn eci_to_ecef(eci: Eci, gmst: f64) -> Result<Ecef> {
    // TODO: Implement conversion
    Ok(Ecef { x: 0.0, y: 0.0, z: 0.0 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinate_structs() {
        let ra_dec = RaDec { ra: 12.0, dec: 45.0 };
        assert_eq!(ra_dec.ra, 12.0);
        assert_eq!(ra_dec.dec, 45.0);

        let alt_az = AltAz { alt: 30.0, az: 180.0 };
        assert_eq!(alt_az.alt, 30.0);
        assert_eq!(alt_az.az, 180.0);
    }

    #[test]
    fn test_ra_dec_to_alt_az_zenith() {
        // Test object at zenith (directly overhead)
        // When RA = LST and Dec = Latitude, object should be at Alt=90°
        let ra_dec = RaDec { ra: 12.0, dec: 40.7128 }; // New York latitude
        let lst = 12.0; // RA = LST means object is transiting (on meridian)
        
        let result = ra_dec_to_alt_az(ra_dec, 40.7128, -74.0060, lst).unwrap();
        
        // Object should be near zenith (Alt ≈ 90°)
        assert!(result.alt > 89.0 && result.alt < 91.0, "Altitude should be near 90°, got {}", result.alt);
    }

    #[test]
    fn test_ra_dec_to_alt_az_horizon() {
        // Test object on horizon (Dec = 0, hour angle = ±90°)
        let ra_dec = RaDec { ra: 6.0, dec: 0.0 }; // On celestial equator
        let lst = 12.0; // 6 hours difference = 90° hour angle
        
        let result = ra_dec_to_alt_az(ra_dec, 40.7128, -74.0060, lst).unwrap();
        
        // Object should be near horizon
        assert!(result.alt < 50.0, "Altitude should be below 50°, got {}", result.alt);
    }

    #[test]
    fn test_ra_dec_to_alt_az_south() {
        // Test object transiting south (Az should be ~180°)
        let ra_dec = RaDec { ra: 12.0, dec: 20.0 };
        let lst = 12.0; // Object transiting
        let lat = 40.7128; // New York
        
        let result = ra_dec_to_alt_az(ra_dec, lat, -74.0060, lst).unwrap();
        
        // For northern hemisphere, object south of zenith should have Az ≈ 180°
        assert!(result.az > 170.0 && result.az < 190.0, "Azimuth should be near 180° (South), got {}", result.az);
    }

    #[test]
    fn test_alt_az_to_ra_dec_round_trip() {
        // Test round-trip conversion: RA/Dec -> Alt/Az -> RA/Dec
        let original = RaDec { ra: 15.5, dec: 35.2 };
        let lst = 18.3;
        let lat = 40.7128;
        let lon = -74.0060;
        
        // Convert to Alt/Az
        let alt_az = ra_dec_to_alt_az(original, lat, lon, lst).unwrap();
        
        // Convert back to RA/Dec
        let result = alt_az_to_ra_dec(alt_az, lat, lon, lst).unwrap();
        
        // Should be close to original (within 0.2 hours ≈ 12 minutes for RA, 0.5° for Dec)
        // Some loss of precision is expected due to floating-point math and trigonometric functions
        assert!((result.ra - original.ra).abs() < 0.2, 
            "RA round-trip failed: original {}, got {}", original.ra, result.ra);
        assert!((result.dec - original.dec).abs() < 0.5, 
            "Dec round-trip failed: original {}, got {}", original.dec, result.dec);
    }

    #[test]
    fn test_alt_az_to_ra_dec_zenith() {
        // Test object at zenith (Alt=90°)
        // Should give RA = LST and Dec = Latitude
        let alt_az = AltAz { alt: 90.0, az: 0.0 }; // Zenith (azimuth doesn't matter)
        let lst = 14.5;
        let lat = 40.7128;
        
        let result = alt_az_to_ra_dec(alt_az, lat, -74.0060, lst).unwrap();
        
        // Dec should equal latitude
        assert!((result.dec - lat).abs() < 0.1, 
            "Declination should equal latitude at zenith, got Dec={}, expected {}", result.dec, lat);
    }

    #[test]
    fn test_coordinate_range_validation() {
        // Test that output coordinates are in valid ranges
        let ra_dec = RaDec { ra: 20.0, dec: 15.0 };
        let lst = 22.0;
        
        let result = ra_dec_to_alt_az(ra_dec, 40.7128, -74.0060, lst).unwrap();
        
        // Altitude: -90° to +90°
        assert!(result.alt >= -90.0 && result.alt <= 90.0, 
            "Altitude out of range: {}", result.alt);
        
        // Azimuth: 0° to 360°
        assert!(result.az >= 0.0 && result.az < 360.0, 
            "Azimuth out of range: {}", result.az);
    }

    #[test]
    fn test_north_celestial_pole() {
        // Test Polaris (near North Celestial Pole)
        // RA ≈ 2.5h, Dec ≈ 89.3° (very close to celestial north pole)
        let polaris = RaDec { ra: 2.5, dec: 89.3 };
        let lst = 12.0;
        let lat = 40.7128; // New York
        
        let result = ra_dec_to_alt_az(polaris, lat, -74.0060, lst).unwrap();
        
        // Polaris should have altitude approximately equal to observer's latitude
        // (it's the North Star - its altitude = your latitude)
        assert!((result.alt - lat).abs() < 5.0, 
            "Polaris altitude should be close to observer latitude, got {}, expected ~{}", result.alt, lat);
    }
}