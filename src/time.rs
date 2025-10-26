//! Time calculations
//! 
//! This module provides functions for astronomical time calculations including:
//! - Julian Date calculations
//! - Sidereal time calculations
//! - Time conversions

use crate::Result;
use chrono::{DateTime, Utc, Datelike, Timelike};

/// Calculate Julian Date from a given date and time
/// 
/// Julian Date is the continuous day count since January 1, 4713 BC at noon UTC.
/// It's the standard time system used in astronomy.
/// 
/// # Arguments
/// * `date_time` - Date and time in UTC
/// 
/// # Returns
/// Julian Date as a floating point number (fractional part represents time of day)
/// 
/// # Example
/// ```
/// use chrono::{DateTime, Utc, NaiveDate, NaiveTime};
/// use cli_astro_calc::time::julian_date;
/// 
/// let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
/// let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
/// let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
/// let jd = julian_date(datetime);
/// // jd should be approximately 2460311.0 for Jan 1, 2024 at noon
/// ```
pub fn julian_date(date_time: DateTime<Utc>) -> f64 {
    let naive = date_time.naive_utc();
    let year = naive.year();
    let month = naive.month();
    let day = naive.day() as f64;
    
    // Get time components as fraction of day
    let hour = naive.hour() as f64;
    let minute = naive.minute() as f64;
    let second = naive.second() as f64;
    let nanosecond = naive.nanosecond() as f64;
    
    // Convert time to fraction of day (0.0 to 1.0)
    // Note: Julian Date starts at noon, so we subtract 12 hours
    let time_fraction = (hour - 12.0 + minute / 60.0 + second / 3600.0 + nanosecond / 3_600_000_000_000.0) / 24.0;
    
    // Julian Date calculation using the standard formula
    let a = (14 - month) / 12;
    let y = year + 4800 - a as i32;
    let m = month + 12 * a - 3;
    
    // Calculate Julian Day Number (at noon)
    let jdn = day + ((153 * m + 2) / 5) as f64 + (365 * y) as f64 + (y / 4) as f64 - (y / 100) as f64 + (y / 400) as f64 - 32045.0;
    
    // Convert to Julian Date (add time fraction from noon)
    jdn + time_fraction
}

/// Calculate Greenwich Mean Sidereal Time (GMST)
/// 
/// Sidereal time is "star time" - Earth's rotation relative to the stars, not the Sun.
/// GMST is the sidereal time at the Greenwich meridian (0° longitude).
/// 
/// # Arguments
/// * `julian_date` - Julian Date
/// 
/// # Returns
/// GMST in hours (0.0 to 24.0)
/// 
/// # Formula
/// GMST = 18.697374558 + 24.06570982441908 * (JD - 2451545.0)
/// Where JD is Julian Date and 2451545.0 is Julian Date for J2000.0 epoch
pub fn greenwich_mean_sidereal_time(julian_date: f64) -> f64 {
    // J2000.0 epoch (January 1, 2000, 12:00:00 UTC)
    const J2000: f64 = 2451545.0;
    
    // Sidereal time rate: 24.06570982441908 hours per Julian day
    const SIDEREAL_RATE: f64 = 24.06570982441908;
    
    // GMST at J2000.0 epoch: 18.697374558 hours
    const GMST_J2000: f64 = 18.697374558;
    
    // Calculate days since J2000.0
    let days_since_j2000 = julian_date - J2000;
    
    // Calculate GMST
    let gmst_hours = GMST_J2000 + SIDEREAL_RATE * days_since_j2000;
    
    // Normalize to 0-24 hour range
    gmst_hours % 24.0
}

/// Calculate Local Sidereal Time (LST)
/// 
/// Local Sidereal Time is the sidereal time at a specific longitude.
/// It's calculated by adding the longitude (converted to time) to GMST.
/// 
/// # Arguments
/// * `gmst` - Greenwich Mean Sidereal Time in hours
/// * `longitude` - Observer longitude in degrees (positive = East, negative = West)
/// 
/// # Returns
/// LST in hours (0.0 to 24.0)
/// 
/// # Formula
/// LST = GMST + longitude / 15.0
/// Where longitude is converted from degrees to hours (15° = 1 hour)
pub fn local_sidereal_time(gmst: f64, longitude: f64) -> f64 {
    // Convert longitude from degrees to hours (15 degrees = 1 hour)
    let longitude_hours = longitude / 15.0;
    
    // Calculate LST
    let lst_hours = gmst + longitude_hours;
    
    // Normalize to 0-24 hour range
    if lst_hours < 0.0 {
        lst_hours + 24.0
    } else if lst_hours >= 24.0 {
        lst_hours - 24.0
    } else {
        lst_hours
    }
}

/// Convert hours to degrees
pub fn hours_to_degrees(hours: f64) -> f64 {
    hours * 15.0
}

/// Convert degrees to hours
pub fn degrees_to_hours(degrees: f64) -> f64 {
    degrees / 15.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::{DateTime, Utc, NaiveDate, NaiveTime};

    #[test]
    fn test_julian_date_j2000() {
        // Test Julian Date for J2000.0 epoch (January 1, 2000, 12:00:00 UTC)
        let date = NaiveDate::from_ymd_opt(2000, 1, 1).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let jd = julian_date(datetime);
        // J2000.0 should be exactly 2451545.0
        assert!((jd - 2451545.0).abs() < 0.001, "J2000.0 Julian Date should be 2451545.0, got {}", jd);
    }

    #[test]
    fn test_julian_date_2024() {
        // Test Julian Date for January 1, 2024, 12:00:00 UTC
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let jd = julian_date(datetime);
        // Should be approximately 2460311.0 (8766 days after J2000.0)
        assert!((jd - 2460311.0).abs() < 0.001, "Jan 1, 2024 Julian Date should be ~2460311.0, got {}", jd);
    }

    #[test]
    fn test_gmst_j2000() {
        // Test GMST at J2000.0 epoch
        let gmst = greenwich_mean_sidereal_time(2451545.0);
        // GMST at J2000.0 should be approximately 18.697374558 hours
        assert!((gmst - 18.697374558).abs() < 0.001, "GMST at J2000.0 should be ~18.697374558, got {}", gmst);
    }

    #[test]
    fn test_lst_calculation() {
        // Test LST calculation for New York (longitude -74.006)
        let gmst = 12.0; // 12:00 GMST
        let longitude = -74.006; // New York longitude
        let lst = local_sidereal_time(gmst, longitude);
        
        // Expected LST = 12.0 + (-74.006/15.0) = 12.0 - 4.9337 = 7.0663
        let expected_lst = 12.0 + (-74.006 / 15.0);
        assert!((lst - expected_lst).abs() < 0.001, "LST calculation incorrect: expected {}, got {}", expected_lst, lst);
    }

    #[test]
    fn test_lst_normalization() {
        // Test LST normalization (should wrap around 24-hour cycle)
        let gmst = 20.0; // 20:00 GMST
        let longitude = 90.0; // 90° East = +6 hours
        let lst = local_sidereal_time(gmst, longitude);
        
        // Expected: 20.0 + 6.0 = 26.0, should normalize to 2.0
        assert!((lst - 2.0).abs() < 0.001, "LST normalization failed: expected 2.0, got {}", lst);
    }

    #[test]
    fn test_hours_degrees_conversion() {
        // Test conversion between hours and degrees
        let hours = 6.0;
        let degrees = hours_to_degrees(hours);
        assert_eq!(degrees, 90.0, "6 hours should equal 90 degrees");
        
        let back_to_hours = degrees_to_hours(degrees);
        assert_eq!(back_to_hours, hours, "Conversion should be reversible");
    }
}
