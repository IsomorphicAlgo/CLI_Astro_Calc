use crate::Result;

#[derive(Debug, Clone, Copy)]
pub struct OrbitalElements {
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub longitude_ascending_node: f64,
    pub argument_periapsis: f64,
    pub mean_anomaly: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct StateVector {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
}

pub fn elements_to_state_vector(elements: OrbitalElements, gm: f64) -> Result<StateVector> {
    let true_anom_rad = mean_to_true_anomaly(elements.mean_anomaly, elements.eccentricity).to_radians();
    let (a, e) = (elements.semi_major_axis, elements.eccentricity);
    let r = a * (1.0 - e * e) / (1.0 + e * true_anom_rad.cos());
    
    let pos_perifocal = [r * true_anom_rad.cos(), r * true_anom_rad.sin(), 0.0];
    let h = (gm * a * (1.0 - e * e)).sqrt();
    let vel_perifocal = [-(gm / h) * true_anom_rad.sin(), (gm / h) * (e + true_anom_rad.cos()), 0.0];
    
    let (incl, raan, argp) = (
        elements.inclination.to_radians(),
        elements.longitude_ascending_node.to_radians(),
        elements.argument_periapsis.to_radians()
    );
    
    let (cr, sr, ci, si, ca, sa) = (raan.cos(), raan.sin(), incl.cos(), incl.sin(), argp.cos(), argp.sin());
    let (r11, r12, r21, r22, r31, r32) = (
        cr * ca - sr * sa * ci,
        -cr * sa - sr * ca * ci,
        sr * ca + cr * sa * ci,
        -sr * sa + cr * ca * ci,
        sa * si,
        ca * si
    );
    
    Ok(StateVector {
        position: [
            r11 * pos_perifocal[0] + r12 * pos_perifocal[1],
            r21 * pos_perifocal[0] + r22 * pos_perifocal[1],
            r31 * pos_perifocal[0] + r32 * pos_perifocal[1],
        ],
        velocity: [
            r11 * vel_perifocal[0] + r12 * vel_perifocal[1],
            r21 * vel_perifocal[0] + r22 * vel_perifocal[1],
            r31 * vel_perifocal[0] + r32 * vel_perifocal[1],
        ],
    })
}

pub fn orbital_period(semi_major_axis: f64, gm: f64) -> f64 {
    use std::f64::consts::PI;
    2.0 * PI * (semi_major_axis.powi(3) / gm).sqrt()
}

pub fn mean_to_true_anomaly(mean_anomaly: f64, eccentricity: f64) -> f64 {
    use std::f64::consts::PI;
    let m_rad = mean_anomaly.to_radians().rem_euclid(2.0 * PI);
    let mut e_anom = if eccentricity < 0.8 { m_rad } else { PI };
    
    for _ in 0..30 {
        let delta = (e_anom - eccentricity * e_anom.sin() - m_rad) / (1.0 - eccentricity * e_anom.cos());
        e_anom -= delta;
        if delta.abs() < 1e-10 { break; }
    }
    
    let true_anom_rad = 2.0 * (((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt() * (e_anom / 2.0).tan()).atan();
    true_anom_rad.to_degrees().rem_euclid(360.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orbital_period_earth() {
        // Earth's orbit around Sun
        // Semi-major axis: ~149,600,000 km
        // Sun's GM: 1.32712440018e11 km³/s²
        let a = 149_600_000.0;
        let gm_sun = 1.32712440018e11;
        
        let period = orbital_period(a, gm_sun);
        
        // Should be approximately 1 year (31,557,600 seconds)
        let one_year_seconds = 365.25 * 24.0 * 3600.0;
        assert!((period - one_year_seconds).abs() < 100_000.0, 
                "Earth's orbital period should be ~1 year, got {} seconds", period);
    }

    #[test]
    fn test_orbital_period_iss() {
        // ISS orbit around Earth
        // Semi-major axis: ~6,780 km (Earth radius ~6,371 km + altitude ~410 km)
        // Earth's GM: 398,600 km³/s²
        let a = 6780.0;
        let gm_earth = 398_600.0;
        
        let period = orbital_period(a, gm_earth);
        
        // ISS orbital period is about 90 minutes (5,400 seconds)
        assert!((period - 5400.0).abs() < 300.0, 
                "ISS orbital period should be ~90 minutes, got {} seconds", period);
    }

    #[test]
    fn test_mean_to_true_anomaly_circular() {
        // For circular orbit (e=0), mean anomaly = true anomaly
        let mean_anom = 45.0;
        let e = 0.0;
        
        let true_anom = mean_to_true_anomaly(mean_anom, e);
        
        assert!((true_anom - mean_anom).abs() < 0.01, 
                "Circular orbit: true anomaly should equal mean anomaly, got {} vs {}", true_anom, mean_anom);
    }

    #[test]
    fn test_mean_to_true_anomaly_eccentric() {
        // Test with eccentric orbit (e=0.5)
        let mean_anom = 90.0;
        let e = 0.5;
        
        let true_anom = mean_to_true_anomaly(mean_anom, e);
        
        // For e=0.5 and M=90°, true anomaly should be > 90° (about 140°)
        // This is because for eccentric orbits, the object spends more time far from periapsis
        assert!(true_anom > 130.0 && true_anom < 150.0, 
                "True anomaly for e=0.5, M=90° should be ~140°, got {}", true_anom);
    }

    #[test]
    fn test_mean_to_true_anomaly_range() {
        // Test that output is always in valid range (0-360)
        let test_cases = vec![0.0, 90.0, 180.0, 270.0, 360.0, 450.0];
        let e = 0.3;
        
        for mean_anom in test_cases {
            let true_anom = mean_to_true_anomaly(mean_anom, e);
            assert!(true_anom >= 0.0 && true_anom < 360.0, 
                    "True anomaly should be in range [0, 360), got {}", true_anom);
        }
    }

    #[test]
    fn test_elements_to_state_vector_circular() {
        // Test circular orbit
        let elements = OrbitalElements {
            semi_major_axis: 7000.0,  // km
            eccentricity: 0.0,
            inclination: 0.0,
            longitude_ascending_node: 0.0,
            argument_periapsis: 0.0,
            mean_anomaly: 0.0,  // At periapsis (which equals apoapsis for circular)
        };
        let gm = 398_600.0;  // Earth's GM
        
        let state = elements_to_state_vector(elements, gm).unwrap();
        
        // Position should be approximately [7000, 0, 0] for mean anomaly = 0
        assert!((state.position[0] - 7000.0).abs() < 1.0, "X position should be ~7000 km");
        assert!(state.position[1].abs() < 1.0, "Y position should be ~0");
        assert!(state.position[2].abs() < 1.0, "Z position should be ~0");
        
        // For circular orbit, velocity magnitude should be sqrt(gm/a)
        let vel_mag = (state.velocity[0].powi(2) + state.velocity[1].powi(2) + state.velocity[2].powi(2)).sqrt();
        let expected_vel = (gm / 7000.0).sqrt();
        assert!((vel_mag - expected_vel).abs() < 0.1, 
                "Velocity magnitude should be {} km/s, got {}", expected_vel, vel_mag);
    }

    #[test]
    fn test_elements_to_state_vector_energy() {
        // Test that specific orbital energy is conserved
        // ε = v²/2 - μ/r = -μ/(2a)
        let elements = OrbitalElements {
            semi_major_axis: 8000.0,
            eccentricity: 0.2,
            inclination: 30.0,
            longitude_ascending_node: 45.0,
            argument_periapsis: 60.0,
            mean_anomaly: 120.0,
        };
        let gm = 398_600.0;
        
        let state = elements_to_state_vector(elements, gm).unwrap();
        
        // Calculate specific energy
        let r_mag = (state.position[0].powi(2) + state.position[1].powi(2) + state.position[2].powi(2)).sqrt();
        let v_mag = (state.velocity[0].powi(2) + state.velocity[1].powi(2) + state.velocity[2].powi(2)).sqrt();
        let energy = v_mag.powi(2) / 2.0 - gm / r_mag;
        let expected_energy = -gm / (2.0 * elements.semi_major_axis);
        
        assert!((energy - expected_energy).abs() < 1.0, 
                "Specific energy should be {} km²/s², got {}", expected_energy, energy);
    }
}
