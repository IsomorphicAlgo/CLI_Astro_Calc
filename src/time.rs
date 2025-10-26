use chrono::{DateTime, Utc, Datelike, Timelike};

pub fn julian_date(date_time: DateTime<Utc>) -> f64 {
    let naive = date_time.naive_utc();
    let (year, month, day) = (naive.year(), naive.month(), naive.day() as f64);
    let (hour, minute, second, nano) = (naive.hour() as f64, naive.minute() as f64, naive.second() as f64, naive.nanosecond() as f64);
    
    let time_frac = (hour - 12.0 + minute / 60.0 + second / 3600.0 + nano / 3_600_000_000_000.0) / 24.0;
    
    let a = (14 - month) / 12;
    let y = year + 4800 - a as i32;
    let m = month + 12 * a - 3;
    
    let jdn = day + ((153 * m + 2) / 5) as f64 + (365 * y) as f64 
              + (y / 4) as f64 - (y / 100) as f64 + (y / 400) as f64 - 32045.0;
    
    jdn + time_frac
}

pub fn greenwich_mean_sidereal_time(julian_date: f64) -> f64 {
    const J2000: f64 = 2451545.0;
    const SIDEREAL_RATE: f64 = 24.06570982441908;
    const GMST_J2000: f64 = 18.697374558;
    
    (GMST_J2000 + SIDEREAL_RATE * (julian_date - J2000)) % 24.0
}

pub fn local_sidereal_time(gmst: f64, longitude: f64) -> f64 {
    (gmst + longitude / 15.0).rem_euclid(24.0)
}

pub fn hours_to_degrees(hours: f64) -> f64 { hours * 15.0 }
pub fn degrees_to_hours(degrees: f64) -> f64 { degrees / 15.0 }

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
