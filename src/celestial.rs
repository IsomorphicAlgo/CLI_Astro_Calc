use crate::Result;
use chrono::{DateTime, Utc};

#[derive(Debug, Clone, Copy)]
pub enum CelestialObject {
    Sun,
    Moon,
}

#[derive(Debug, Clone, Copy)]
pub struct ObserverLocation {
    pub latitude: f64,
    pub longitude: f64,
    pub elevation: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct RiseSetTimes {
    pub rise: Option<DateTime<Utc>>,
    pub set: Option<DateTime<Utc>>,
}

pub fn calculate_rise_set_times(
    object: CelestialObject,
    location: ObserverLocation,
    date: DateTime<Utc>,
) -> Result<RiseSetTimes> {
    match object {
        CelestialObject::Sun => calculate_solar_rise_set(location, date),
        CelestialObject::Moon => calculate_lunar_rise_set(location, date),
    }
}

fn calculate_solar_rise_set(location: ObserverLocation, date: DateTime<Utc>) -> Result<RiseSetTimes> {
    use crate::time::{julian_date, local_sidereal_time, greenwich_mean_sidereal_time};
    use chrono::NaiveTime;
    
    let naive_date = date.naive_utc().date();
    let midnight = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap()), chrono::Utc);
    let noon = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(12, 0, 0).unwrap()), chrono::Utc);
    
    let solar_pos = calculate_solar_position(noon)?;
    let dec_rad = solar_pos.dec.to_radians();
    let lat_rad = location.latitude.to_radians();
    let altitude_correction = -0.833_f64.to_radians();
    
    let cos_hour_angle = (altitude_correction.sin() - lat_rad.sin() * dec_rad.sin()) 
                       / (lat_rad.cos() * dec_rad.cos());
    
    if cos_hour_angle.abs() > 1.0 {
        return Ok(RiseSetTimes { rise: None, set: None });
    }
    
    let hour_angle_deg = cos_hour_angle.acos().to_degrees();
    let lst_midnight = local_sidereal_time(
        greenwich_mean_sidereal_time(julian_date(midnight)), location.longitude);
    let transit_time = (solar_pos.ra - lst_midnight).rem_euclid(24.0);
    
    let rise_time_hours = (transit_time - hour_angle_deg / 15.0).rem_euclid(24.0);
    let set_time_hours = (transit_time + hour_angle_deg / 15.0).rem_euclid(24.0);
    
    Ok(RiseSetTimes {
        rise: Some(midnight + chrono::Duration::seconds((rise_time_hours * 3600.0) as i64)),
        set: Some(midnight + chrono::Duration::seconds((set_time_hours * 3600.0) as i64)),
    })
}

fn calculate_lunar_rise_set(location: ObserverLocation, date: DateTime<Utc>) -> Result<RiseSetTimes> {
    use crate::time::{julian_date, local_sidereal_time, greenwich_mean_sidereal_time};
    use chrono::NaiveTime;
    
    let naive_date = date.naive_utc().date();
    let midnight = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(0, 0, 0).unwrap()), chrono::Utc);
    let noon = DateTime::from_naive_utc_and_offset(
        naive_date.and_time(NaiveTime::from_hms_opt(12, 0, 0).unwrap()), chrono::Utc);
    
    let lunar_pos = calculate_lunar_position(noon)?;
    let dec_rad = lunar_pos.dec.to_radians();
    let lat_rad = location.latitude.to_radians();
    let altitude_correction = -0.583_f64.to_radians();
    
    let cos_hour_angle = (altitude_correction.sin() - lat_rad.sin() * dec_rad.sin()) 
                       / (lat_rad.cos() * dec_rad.cos());
    
    if cos_hour_angle.abs() > 1.0 {
        return Ok(RiseSetTimes { rise: None, set: None });
    }
    
    let hour_angle_deg = cos_hour_angle.acos().to_degrees();
    let lst_midnight = local_sidereal_time(
        greenwich_mean_sidereal_time(julian_date(midnight)), location.longitude);
    let transit_time = (lunar_pos.ra - lst_midnight).rem_euclid(24.0);
    
    let rise_time_hours = (transit_time - hour_angle_deg / 15.0).rem_euclid(24.0);
    let set_time_hours = (transit_time + hour_angle_deg / 15.0).rem_euclid(24.0);
    
    Ok(RiseSetTimes {
        rise: Some(midnight + chrono::Duration::seconds((rise_time_hours * 3600.0) as i64)),
        set: Some(midnight + chrono::Duration::seconds((set_time_hours * 3600.0) as i64)),
    })
}

pub fn calculate_position(object: CelestialObject, date: DateTime<Utc>) -> Result<crate::coordinates::RaDec> {
    match object {
        CelestialObject::Sun => calculate_solar_position(date),
        CelestialObject::Moon => calculate_lunar_position(date),
    }
}

fn calculate_solar_position(date: DateTime<Utc>) -> Result<crate::coordinates::RaDec> {
    use crate::time::julian_date;
    
    let jd = julian_date(date);
    const J2000: f64 = 2451545.0;
    let d = jd - J2000;
    
    let mean_anomaly = 357.5291 + 0.98560028 * d;
    let m_rad = mean_anomaly.to_radians();
    let eoc = 1.9148 * m_rad.sin() + 0.0200 * (2.0 * m_rad).sin() + 0.0003 * (3.0 * m_rad).sin();
    
    let ecliptic_lon = mean_anomaly + eoc + 180.0 + 102.9372;
    let obliquity = 23.4393 - 0.0000004 * d;
    
    let lambda = ecliptic_lon.to_radians();
    let epsilon = obliquity.to_radians();
    
    let ra_rad = (lambda.sin() * epsilon.cos()).atan2(lambda.cos());
    let dec_rad = (lambda.sin() * epsilon.sin()).asin();
    
    Ok(crate::coordinates::RaDec {
        ra: (ra_rad.to_degrees() / 15.0).rem_euclid(24.0),
        dec: dec_rad.to_degrees(),
    })
}

fn calculate_lunar_position(date: DateTime<Utc>) -> Result<crate::coordinates::RaDec> {
    use crate::time::julian_date;
    
    let jd = julian_date(date);
    const J2000: f64 = 2451545.0;
    let t = (jd - J2000) / 36525.0;
    
    let l_prime = 218.3164477 + 481267.88123421 * t - 0.0015786 * t * t + t.powi(3) / 538841.0 - t.powi(4) / 65194000.0;
    let d = 297.8501921 + 445267.1114034 * t - 0.0018819 * t * t + t.powi(3) / 545868.0 - t.powi(4) / 113065000.0;
    let m = 357.5291092 + 35999.0502909 * t - 0.0001536 * t * t + t.powi(3) / 24490000.0;
    let m_prime = 134.9633964 + 477198.8675055 * t + 0.0087414 * t * t + t.powi(3) / 69699.0 - t.powi(4) / 14712000.0;
    let f = 93.2720950 + 483202.0175233 * t - 0.0036539 * t * t - t.powi(3) / 3526000.0 + t.powi(4) / 863310000.0;
    
    let (dr, mr, mpr, fr) = (d.to_radians(), m.to_radians(), m_prime.to_radians(), f.to_radians());
    
    let sigma_l = 6288774.0 * mpr.sin() + 1274027.0 * (2.0 * dr - mpr).sin() + 658314.0 * (2.0 * dr).sin()
                + 213618.0 * (2.0 * mpr).sin() - 185116.0 * mr.sin() - 114332.0 * (2.0 * fr).sin()
                + 58793.0 * (2.0 * dr - 2.0 * mpr).sin() + 57066.0 * (2.0 * dr - mr - mpr).sin()
                + 53322.0 * (2.0 * dr + mpr).sin() + 45758.0 * (2.0 * dr - mr).sin();
    
    let sigma_b = 5128122.0 * fr.sin() + 280602.0 * (mpr + fr).sin() + 277693.0 * (mpr - fr).sin()
                + 173237.0 * (2.0 * dr - fr).sin() + 55413.0 * (2.0 * dr - mpr + fr).sin()
                + 46271.0 * (2.0 * dr - mpr - fr).sin() + 32573.0 * (2.0 * dr + fr).sin()
                + 17198.0 * (2.0 * mpr + fr).sin();
    
    let lambda = (l_prime + sigma_l / 1000000.0).rem_euclid(360.0);
    let beta = (sigma_b / 1000000.0).rem_euclid(360.0);
    let epsilon = 23.439291 - 0.0130042 * t - 0.00000016 * t * t + 0.000000504 * t.powi(3);
    
    let (lambda_rad, beta_rad, epsilon_rad) = (lambda.to_radians(), beta.to_radians(), epsilon.to_radians());
    
    let ra_rad = (lambda_rad.sin() * epsilon_rad.cos() - beta_rad.tan() * epsilon_rad.sin()).atan2(lambda_rad.cos());
    let dec_rad = (beta_rad.sin() * epsilon_rad.cos() + beta_rad.cos() * epsilon_rad.sin() * lambda_rad.sin()).asin();
    
    Ok(crate::coordinates::RaDec {
        ra: (ra_rad.to_degrees() / 15.0).rem_euclid(24.0),
        dec: dec_rad.to_degrees(),
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

    #[test]
    fn test_lunar_position_basic() {
        // Test lunar position calculation at a known date
        let date = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let time = NaiveTime::from_hms_opt(0, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_lunar_position(datetime).unwrap();
        
        // Verify RA is in valid range (0-24 hours)
        assert!(position.ra >= 0.0 && position.ra < 24.0, "Lunar RA should be in range 0-24 hours, got {}", position.ra);
        
        // Verify Dec is in valid range (-90° to +90°)
        assert!(position.dec >= -90.0 && position.dec <= 90.0, "Lunar Dec should be in range -90° to +90°, got {}", position.dec);
    }

    #[test]
    fn test_lunar_position_declination_range() {
        // Moon's declination should be within ±28.5° (maximum lunar declination)
        let date = NaiveDate::from_ymd_opt(2024, 6, 15).unwrap();
        let time = NaiveTime::from_hms_opt(12, 0, 0).unwrap();
        let datetime = DateTime::from_naive_utc_and_offset(date.and_time(time), Utc);
        
        let position = calculate_lunar_position(datetime).unwrap();
        
        // Moon's orbit is inclined about 5° to ecliptic, which is inclined 23.5° to equator
        // So max declination is about ±28.5°
        assert!(position.dec.abs() <= 29.0, "Lunar declination should be within ±29°, got {}", position.dec);
    }

    #[test]
    fn test_lunar_position_changes_over_month() {
        // Moon should move significantly over a month (orbital period ~27.3 days)
        let date1 = NaiveDate::from_ymd_opt(2024, 1, 1).unwrap();
        let time1 = NaiveTime::from_hms_opt(0, 0, 0).unwrap();
        let datetime1 = DateTime::from_naive_utc_and_offset(date1.and_time(time1), Utc);
        
        let date2 = NaiveDate::from_ymd_opt(2024, 1, 28).unwrap();
        let time2 = NaiveTime::from_hms_opt(0, 0, 0).unwrap();
        let datetime2 = DateTime::from_naive_utc_and_offset(date2.and_time(time2), Utc);
        
        let pos1 = calculate_lunar_position(datetime1).unwrap();
        let pos2 = calculate_lunar_position(datetime2).unwrap();
        
        // Verify positions are different (Moon moves ~13° per day)
        // Over 27 days, should complete almost full orbit
        let ra_diff = (pos2.ra - pos1.ra).abs();
        
        // Just verify the position has changed (any non-zero difference is good)
        // The Moon's actual motion is complex due to perturbations
        assert!(ra_diff > 0.01, 
                "Moon position should change over a month, RA1: {}, RA2: {}", pos1.ra, pos2.ra);
    }
}
