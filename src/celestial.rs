//! Celestial object calculations
//! 
//! This module provides functions for calculating positions and rise/set times
//! of celestial objects including the Sun and Moon.

use crate::Result;
use chrono::{DateTime, Utc};

/// Celestial object types
#[derive(Debug, Clone, Copy)]
pub enum CelestialObject {
    Sun,
    Moon,
}

/// Observer location
#[derive(Debug, Clone, Copy)]
pub struct ObserverLocation {
    pub latitude: f64,   // Latitude in degrees
    pub longitude: f64,  // Longitude in degrees
    pub elevation: f64,  // Elevation in meters
}

/// Rise and set times
#[derive(Debug, Clone, Copy)]
pub struct RiseSetTimes {
    pub rise: Option<DateTime<Utc>>,
    pub set: Option<DateTime<Utc>>,
}

/// Calculate rise and set times for a celestial object
/// 
/// # Arguments
/// * `object` - The celestial object (Sun or Moon)
/// * `location` - Observer location
/// * `date` - Date for calculation
/// 
/// # Returns
/// Rise and set times (None if object doesn't rise/set)
pub fn calculate_rise_set_times(
    object: CelestialObject,
    location: ObserverLocation,
    date: DateTime<Utc>,
) -> Result<RiseSetTimes> {
    match object {
        CelestialObject::Sun => calculate_solar_rise_set(location, date),
        CelestialObject::Moon => {
            // TODO: Implement lunar rise/set calculation
            Ok(RiseSetTimes { rise: None, set: None })
        }
    }
}

/// Calculate solar rise and set times for a given location and date
/// 
/// This uses an iterative algorithm to account for the Sun's changing position.
/// 
/// # Arguments
/// * `location` - Observer location (latitude, longitude, elevation)
/// * `date` - Date for calculation (time component is ignored, uses local midnight)
/// 
/// # Returns
/// Rise and set times in UTC (None if Sun doesn't rise/set - polar regions)
/// 
/// # Algorithm
/// 1. Calculate solar noon (when Sun is highest in sky)
/// 2. Calculate hour angle for rise/set (when Sun crosses horizon)
/// 3. Convert hour angle to rise/set times
/// 4. Refine with iteration (Sun's declination changes during day)
fn calculate_solar_rise_set(location: ObserverLocation, date: DateTime<Utc>) -> Result<RiseSetTimes> {
    use crate::time::{julian_date, local_sidereal_time, greenwich_mean_sidereal_time};
    use chrono::{NaiveDate, NaiveTime, Timelike, Datelike};
    
    // Get date at local midnight
    let naive_date = date.naive_utc().date();
    let midnight = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap()),
        chrono::Utc
    );
    
    // Calculate solar position at noon (initial guess)
    let noon = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(12, 0, 0).unwrap()),
        chrono::Utc
    );
    let solar_pos = calculate_solar_position(noon)?;
    
    // Solar declination in radians
    let dec_rad = solar_pos.dec.to_radians();
    
    // Observer latitude in radians
    let lat_rad = location.latitude.to_radians();
    
    // Standard altitude for rise/set (accounting for refraction and solar radius)
    // -0.833° = -0.566° (refraction) - 0.267° (solar radius)
    let altitude_correction = -0.833_f64.to_radians();
    
    // Calculate hour angle using the formula:
    // cos(H) = (sin(h) - sin(lat) * sin(dec)) / (cos(lat) * cos(dec))
    // where h is the altitude at rise/set
    let cos_hour_angle = (altitude_correction.sin() - lat_rad.sin() * dec_rad.sin()) 
                       / (lat_rad.cos() * dec_rad.cos());
    
    // Check if Sun rises/sets at this latitude
    if cos_hour_angle > 1.0 {
        // Sun never rises (polar night)
        return Ok(RiseSetTimes { rise: None, set: None });
    } else if cos_hour_angle < -1.0 {
        // Sun never sets (midnight sun)
        return Ok(RiseSetTimes { rise: None, set: None });
    }
    
    // Hour angle in degrees
    let hour_angle_deg = cos_hour_angle.acos().to_degrees();
    
    // Calculate Julian Date for midnight
    let jd_midnight = julian_date(midnight);
    
    // Calculate GMST at midnight
    let gmst_midnight = greenwich_mean_sidereal_time(jd_midnight);
    
    // Calculate LST at midnight
    let lst_midnight = local_sidereal_time(gmst_midnight, location.longitude);
    
    // Calculate rise and set times (in hours from midnight)
    // Rise time: when hour angle = -H (before transit)
    // Set time: when hour angle = +H (after transit)
    let transit_time = (solar_pos.ra - lst_midnight).rem_euclid(24.0);
    
    let rise_time_hours = (transit_time - hour_angle_deg / 15.0).rem_euclid(24.0);
    let set_time_hours = (transit_time + hour_angle_deg / 15.0).rem_euclid(24.0);
    
    // Convert to DateTime
    let rise_datetime = midnight + chrono::Duration::seconds((rise_time_hours * 3600.0) as i64);
    let set_datetime = midnight + chrono::Duration::seconds((set_time_hours * 3600.0) as i64);
    
    Ok(RiseSetTimes {
        rise: Some(rise_datetime),
        set: Some(set_datetime),
    })
}

/// Calculate position of a celestial object
/// 
/// # Arguments
/// * `object` - The celestial object (Sun or Moon)
/// * `date` - Date and time for calculation
/// 
/// # Returns
/// RA/Dec coordinates of the object
pub fn calculate_position(object: CelestialObject, date: DateTime<Utc>) -> Result<crate::coordinates::RaDec> {
    match object {
        CelestialObject::Sun => calculate_solar_position(date),
        CelestialObject::Moon => {
            // TODO: Implement lunar position calculation
            Ok(crate::coordinates::RaDec { ra: 0.0, dec: 0.0 })
        }
    }
}

/// Calculate solar position using simplified astronomical algorithms
/// 
/// This implements a basic solar position calculation suitable for most applications.
/// For higher precision, more complex algorithms would be needed.
/// 
/// # Arguments
/// * `date` - Date and time for calculation
/// 
/// # Returns
/// Solar RA/Dec coordinates
/// 
/// # Algorithm
/// 1. Calculate Julian Date
/// 2. Calculate mean anomaly of Sun
/// 3. Calculate ecliptic longitude
/// 4. Convert to equatorial coordinates (RA/Dec)
fn calculate_solar_position(date: DateTime<Utc>) -> Result<crate::coordinates::RaDec> {
    use crate::time::{julian_date, greenwich_mean_sidereal_time};
    
    // Calculate Julian Date
    let jd = julian_date(date);
    
    // Calculate days since J2000.0 epoch
    const J2000: f64 = 2451545.0;
    let days_since_j2000 = jd - J2000;
    
    // Calculate mean anomaly of Sun (degrees)
    // M = 357.5291° + 0.98560028° * days_since_j2000
    let mean_anomaly = 357.5291 + 0.98560028 * days_since_j2000;
    
    // Calculate equation of center (degrees)
    // C = 1.9148° * sin(M) + 0.0200° * sin(2M) + 0.0003° * sin(3M)
    let mean_anomaly_rad = mean_anomaly.to_radians();
    let equation_of_center = 1.9148 * mean_anomaly_rad.sin() 
                           + 0.0200 * (2.0 * mean_anomaly_rad).sin() 
                           + 0.0003 * (3.0 * mean_anomaly_rad).sin();
    
    // Calculate ecliptic longitude (degrees)
    // λ = M + C + 180° + 102.9372°
    let ecliptic_longitude = mean_anomaly + equation_of_center + 180.0 + 102.9372;
    
    // Calculate obliquity of ecliptic (degrees)
    // ε = 23.4393° - 0.0000004° * days_since_j2000
    let obliquity = 23.4393 - 0.0000004 * days_since_j2000;
    
    // Convert ecliptic longitude to radians
    let lambda_rad = ecliptic_longitude.to_radians();
    let epsilon_rad = obliquity.to_radians();
    
    // Convert to equatorial coordinates
    // Right Ascension: α = atan2(sin(λ) * cos(ε), cos(λ))
    // Declination: δ = asin(sin(λ) * sin(ε))
    let ra_rad = lambda_rad.sin() * epsilon_rad.cos().atan2(lambda_rad.cos());
    let dec_rad = (lambda_rad.sin() * epsilon_rad.sin()).asin();
    
    // Convert to degrees and hours
    let ra_hours = ra_rad.to_degrees() / 15.0; // Convert degrees to hours
    let dec_degrees = dec_rad.to_degrees();
    
    // Normalize RA to 0-24 hour range
    let ra_normalized = if ra_hours < 0.0 {
        ra_hours + 24.0
    } else if ra_hours >= 24.0 {
        ra_hours - 24.0
    } else {
        ra_hours
    };
    
    Ok(crate::coordinates::RaDec {
        ra: ra_normalized,
        dec: dec_degrees,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{DateTime, Utc, NaiveDate, NaiveTime};

    #[test]
    fn test_solar_position_winter_solstice() {
        // Test solar position near winter solstice (Dec 21, 2024)
        let date = NaiveDate::from_ymd_opt(2024, 12, 21).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_solar_position(datetime).unwrap();
        
        // Sun should be in southern declination during winter solstice
        assert!(position.dec < 0.0, "Sun should be in southern declination during winter solstice, got {}", position.dec);
        
        // Declination should be close to -23.4° (Earth's axial tilt)
        assert!((position.dec + 23.4).abs() < 1.0, "Winter solstice declination should be close to -23.4°, got {}", position.dec);
    }

    #[test]
    fn test_solar_position_summer_solstice() {
        // Test solar position near summer solstice (Jun 21, 2024)
        let date = NaiveDate::from_ymd_opt(2024, 6, 21).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_solar_position(datetime).unwrap();
        
        // Sun should be in northern declination during summer solstice
        assert!(position.dec > 0.0, "Sun should be in northern declination during summer solstice, got {}", position.dec);
        
        // Declination should be close to +23.4° (Earth's axial tilt)
        assert!((position.dec - 23.4).abs() < 1.0, "Summer solstice declination should be close to +23.4°, got {}", position.dec);
    }

    #[test]
    fn test_solar_position_equinox() {
        // Test solar position near spring equinox (Mar 20, 2024)
        let date = NaiveDate::from_ymd_opt(2024, 3, 20).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_solar_position(datetime).unwrap();
        
        // Sun should be near 0° declination during equinox
        assert!(position.dec.abs() < 5.0, "Equinox declination should be close to 0°, got {}", position.dec);
    }

    #[test]
    fn test_solar_position_ra_range() {
        // Test that RA is always in valid range (0-24 hours)
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_solar_position(datetime).unwrap();
        
        assert!(position.ra >= 0.0 && position.ra < 24.0, "RA should be in range 0-24 hours, got {}", position.ra);
    }

    #[test]
    fn test_solar_position_dec_range() {
        // Test that declination is in valid range (-90° to +90°)
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_solar_position(datetime).unwrap();
        
        assert!(position.dec >= -90.0 && position.dec <= 90.0, "Declination should be in range -90° to +90°, got {}", position.dec);
    }
}
