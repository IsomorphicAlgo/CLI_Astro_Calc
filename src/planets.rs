use crate::Result;
use crate::coordinates::RaDec;

/// Represents the major planets in the solar system.
/// 
/// Planets are ordered by distance from the Sun (excluding Earth, which is the observer).
/// Earth is included for completeness but typically not used for position calculations
/// from Earth's perspective (use geocentric calculations instead).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Planet {
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
}

impl Planet {
    /// Returns the name of the planet as a string.
    pub fn name(&self) -> &'static str {
        match self {
            Planet::Mercury => "Mercury",
            Planet::Venus => "Venus",
            Planet::Earth => "Earth",
            Planet::Mars => "Mars",
            Planet::Jupiter => "Jupiter",
            Planet::Saturn => "Saturn",
            Planet::Uranus => "Uranus",
            Planet::Neptune => "Neptune",
        }
    }

    /// Parses a planet name from a string (case-insensitive).
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "mercury" => Some(Planet::Mercury),
            "venus" => Some(Planet::Venus),
            "earth" => Some(Planet::Earth),
            "mars" => Some(Planet::Mars),
            "jupiter" => Some(Planet::Jupiter),
            "saturn" => Some(Planet::Saturn),
            "uranus" => Some(Planet::Uranus),
            "neptune" => Some(Planet::Neptune),
            _ => None,
        }
    }
}

/// VSOP87 coefficient for a single term in the series.
/// 
/// Each term in VSOP87 has the form: A × cos(B + C × t)
/// where:
/// - A: Amplitude (coefficient)
/// - B: Phase (constant term in radians)
/// - C: Frequency (coefficient for time variable)
/// - t: Time in Julian centuries from J2000.0
#[derive(Debug, Clone, Copy)]
pub struct Vsop87Term {
    /// Amplitude A
    pub amplitude: f64,
    /// Phase B (in radians)
    pub phase: f64,
    /// Frequency C
    pub frequency: f64,
}

/// VSOP87 series for a single variable (L, B, or R).
/// 
/// VSOP87 represents each coordinate as a sum of series:
/// - L0, L1, L2, L3, L4, L5 for longitude
/// - B0, B1, B2, B3, B4 for latitude
/// - R0, R1, R2, R3, R4 for radius
/// 
/// Each series is a sum of terms: Σ(A × cos(B + C × t))
/// The final value is: (L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8
#[derive(Debug, Clone)]
pub struct Vsop87Series {
    /// Series L0, B0, or R0 (constant term)
    pub series_0: Vec<Vsop87Term>,
    /// Series L1, B1, or R1 (linear term)
    pub series_1: Vec<Vsop87Term>,
    /// Series L2, B2, or R2 (quadratic term)
    pub series_2: Vec<Vsop87Term>,
    /// Series L3, B3, or R3 (cubic term)
    pub series_3: Vec<Vsop87Term>,
    /// Series L4, B4, or R4 (quartic term) - only for L
    pub series_4: Option<Vec<Vsop87Term>>,
    /// Series L5 (quintic term) - only for L
    pub series_5: Option<Vec<Vsop87Term>>,
}

/// Complete VSOP87 data for a planet.
/// 
/// Contains three VSOP87Series: one for longitude (L), one for latitude (B),
/// and one for radius (R).
#[derive(Debug, Clone)]
pub struct PlanetVsop87Data {
    /// Longitude series (L0, L1, L2, L3, L4, L5)
    pub longitude: Vsop87Series,
    /// Latitude series (B0, B1, B2, B3, B4)
    pub latitude: Vsop87Series,
    /// Radius series (R0, R1, R2, R3, R4) in AU
    pub radius: Vsop87Series,
}

/// Heliocentric ecliptic coordinates from VSOP87.
/// 
/// These are the raw output from VSOP87 calculations before conversion
/// to equatorial coordinates (RA/Dec).
#[derive(Debug, Clone, Copy)]
pub struct HeliocentricEcliptic {
    /// Ecliptic longitude (L) in radians
    pub longitude: f64,
    /// Ecliptic latitude (B) in radians
    pub latitude: f64,
    /// Radius vector (R) in Astronomical Units (AU)
    pub radius: f64,
}

/// Calculates the position of a planet using VSOP87 theory.
/// 
/// This function computes the heliocentric ecliptic coordinates (L, B, R)
/// and converts them to geocentric equatorial coordinates (RA, Dec).
/// 
/// # Arguments
/// * `planet` - The planet to calculate
/// * `julian_date` - Julian Date for the calculation
/// 
/// # Returns
/// Right Ascension and Declination in equatorial coordinates
/// 
/// # Errors
/// Returns an error if:
/// - Planet data is not available
/// - Julian Date is invalid
/// - Calculation fails
/// 
/// # Accuracy
/// For a simplified VSOP87 implementation (truncated series):
/// - Inner planets (Mercury, Venus, Mars): ~1 arcminute
/// - Outer planets (Jupiter, Saturn): ~1-2 arcminutes
/// - Distant planets (Uranus, Neptune): ~2-5 arcminutes
/// 
/// Full VSOP87 accuracy is typically better than 1 arcsecond for all planets.
/// 
/// # Example
/// ```rust
/// use cli_astro_calc::planets::{Planet, calculate_planet_position};
/// 
/// // Calculate Jupiter's position on January 1, 2000
/// let jd = 2451545.0; // J2000.0
/// let position = calculate_planet_position(Planet::Jupiter, jd)?;
/// println!("Jupiter: RA={:.2}h, Dec={:.2}°", position.ra, position.dec);
/// ```
pub fn calculate_planet_position(planet: Planet, _julian_date: f64) -> Result<RaDec> {
    // Get VSOP87 data for planet
    let _vsop87_data = get_planet_vsop87_data(planet)
        .ok_or_else(|| crate::error::AstroError::InvalidCoordinate(
            format!("VSOP87 data not available for {}", planet.name())
        ))?;
    
    // TODO: Implement in Step B3
    // 2. Calculate time in Julian centuries from J2000.0
    // 3. Evaluate VSOP87 series for L, B, R using vsop87_data
    // 4. Convert heliocentric ecliptic to geocentric equatorial
    // 5. Return RA/Dec
    
    Err(crate::error::AstroError::InvalidCoordinate(
        format!("Planet position calculation not yet fully implemented for {}. VSOP87 data structure ready.", planet.name())
    ))
}

/// Evaluates a VSOP87 series for a given time.
/// 
/// # Arguments
/// * `series` - The VSOP87 series to evaluate
/// * `t` - Time in Julian centuries from J2000.0
/// 
/// # Returns
/// The evaluated series value (in radians for L/B, AU for R)
/// 
/// # Formula
/// result = (L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8
/// where each Ln = Σ(A × cos(B + C × t))
fn evaluate_vsop87_series(_series: &Vsop87Series, _t: f64) -> f64 {
    // TODO: Implement in Step B3
    // Evaluate each series (0-5) and combine with time powers
    0.0
}

/// Evaluates a single VSOP87 term: A × cos(B + C × t)
/// 
/// # Arguments
/// * `term` - The VSOP87 term
/// * `t` - Time in Julian centuries from J2000.0
/// 
/// # Returns
/// The evaluated term value
fn evaluate_vsop87_term(term: &Vsop87Term, t: f64) -> f64 {
    term.amplitude * (term.phase + term.frequency * t).cos()
}

/// Converts heliocentric ecliptic coordinates to geocentric equatorial coordinates.
/// 
/// This conversion accounts for:
/// 1. Earth's position in the solar system
/// 2. Ecliptic to equatorial coordinate transformation
/// 3. Light-time correction (optional, for high precision)
/// 
/// # Arguments
/// * `heliocentric` - Heliocentric ecliptic coordinates (L, B, R)
/// * `earth_position` - Earth's heliocentric position (for geocentric conversion)
/// * `julian_date` - Julian Date for obliquity calculation
/// 
/// # Returns
/// Geocentric equatorial coordinates (RA, Dec)
/// 
/// # Note
/// For initial implementation, light-time correction may be omitted.
/// This reduces accuracy by ~0.01 arcseconds but simplifies the calculation.
fn heliocentric_to_geocentric(
    _heliocentric: HeliocentricEcliptic,
    _earth_position: HeliocentricEcliptic,
    _julian_date: f64,
) -> Result<RaDec> {
    // TODO: Implement in Step B4
    // 1. Calculate geocentric position vector
    // 2. Convert ecliptic to equatorial using obliquity
    // 3. Convert to RA/Dec
    // 4. Return result
    
    Err(crate::error::AstroError::InvalidCoordinate(
        "Geocentric conversion not yet implemented".to_string()
    ))
}

/// Gets VSOP87 data for a planet.
/// 
/// # Arguments
/// * `planet` - The planet to get data for
/// 
/// # Returns
/// VSOP87 coefficient data for the planet
/// 
/// # Note
/// This implementation uses simplified/truncated VSOP87 data.
/// Full VSOP87 data contains thousands of terms per planet.
/// The truncated version uses only the most significant terms for each series,
/// providing ~1-5 arcminute accuracy (vs <1 arcsecond for full VSOP87).
/// 
/// # Data Source
/// Coefficients are based on VSOP87 theory (Bretagnon & Francou, 1987).
/// Truncated versions use only terms with significant amplitudes.
/// For full precision, download complete VSOP87 data from IMCCE:
/// https://ftp.imcce.fr/pub/ephem/planets/vsop87/
pub(crate) fn get_planet_vsop87_data(planet: Planet) -> Option<PlanetVsop87Data> {
    match planet {
        Planet::Mercury => Some(get_mercury_vsop87_data()),
        Planet::Venus => Some(get_venus_vsop87_data()),
        Planet::Earth => Some(get_earth_vsop87_data()),
        Planet::Mars => Some(get_mars_vsop87_data()),
        Planet::Jupiter => Some(get_jupiter_vsop87_data()),
        Planet::Saturn => Some(get_saturn_vsop87_data()),
        Planet::Uranus => Some(get_uranus_vsop87_data()),
        Planet::Neptune => Some(get_neptune_vsop87_data()),
    }
}

/// Gets simplified VSOP87 data for Mercury.
/// 
/// This is a truncated version using only the most significant terms.
/// For demonstration purposes, includes a minimal set of coefficients.
fn get_mercury_vsop87_data() -> PlanetVsop87Data {
    // Simplified VSOP87 data for Mercury
    // These are truncated coefficients - only most significant terms included
    // Full VSOP87 for Mercury has hundreds of terms per series
    
    // Longitude (L) series - truncated to largest terms
    let longitude = Vsop87Series {
        // L0: Constant term series (largest terms only)
        series_0: vec![
            Vsop87Term { amplitude: 4.40250710144, phase: 0.0, frequency: 0.0 },
            Vsop87Term { amplitude: 0.40989414976, phase: 1.48302034194, frequency: 26087.90314157420 },
            Vsop87Term { amplitude: 0.05046294199, phase: 4.47785489540, frequency: 52175.80628314840 },
            Vsop87Term { amplitude: 0.00855346843, phase: 1.16520322359, frequency: 78263.70942462180 },
            Vsop87Term { amplitude: 0.00165506162, phase: 4.11969133181, frequency: 104351.61256629680 },
        ],
        // L1: Linear term series
        series_1: vec![
            Vsop87Term { amplitude: 26087.90314157420, phase: 0.0, frequency: 0.0 },
            Vsop87Term { amplitude: 0.01126007832, phase: 6.21703970996, frequency: 26087.90314157420 },
            Vsop87Term { amplitude: 0.00303471395, phase: 3.05524609620, frequency: 52175.80628314840 },
        ],
        // L2: Quadratic term series
        series_2: vec![
            Vsop87Term { amplitude: 0.00053049845, phase: 0.0, frequency: 0.0 },
            Vsop87Term { amplitude: 0.00016903658, phase: 4.69072300649, frequency: 26087.90314157420 },
        ],
        // L3: Cubic term series
        series_3: vec![
            Vsop87Term { amplitude: 0.00000169496, phase: 3.20221586859, frequency: 26087.90314157420 },
        ],
        // L4: Quartic term series (optional, only for L)
        series_4: Some(vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 }, // Minimal term
        ]),
        // L5: Quintic term series (optional, only for L)
        series_5: Some(vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 }, // Minimal term
        ]),
    };
    
    // Latitude (B) series - truncated
    let latitude = Vsop87Series {
        series_0: vec![
            Vsop87Term { amplitude: 0.11737528962, phase: 1.98357498767, frequency: 26087.90314157420 },
            Vsop87Term { amplitude: 0.02388076996, phase: 5.03738959686, frequency: 52175.80628314840 },
            Vsop87Term { amplitude: 0.01222839532, phase: 3.14159265359, frequency: 0.0 },
        ],
        series_1: vec![
            Vsop87Term { amplitude: 0.00397535498, phase: 4.93750888835, frequency: 26087.90314157420 },
        ],
        series_2: vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 },
        ],
        series_3: vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 },
        ],
        series_4: None, // B4 not used
        series_5: None, // B5 not used
    };
    
    // Radius (R) series - truncated
    let radius = Vsop87Series {
        series_0: vec![
            Vsop87Term { amplitude: 0.39528271652, phase: 0.0, frequency: 0.0 },
            Vsop87Term { amplitude: 0.07834131717, phase: 6.19233722599, frequency: 26087.90314157420 },
            Vsop87Term { amplitude: 0.00795532757, phase: 2.95989680096, frequency: 52175.80628314840 },
        ],
        series_1: vec![
            Vsop87Term { amplitude: 0.00217347739, phase: 4.65617158663, frequency: 26087.90314157420 },
        ],
        series_2: vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 },
        ],
        series_3: vec![
            Vsop87Term { amplitude: 0.00000000000, phase: 0.0, frequency: 0.0 },
        ],
        series_4: None, // R4 not used
        series_5: None, // R5 not used
    };
    
    PlanetVsop87Data {
        longitude,
        latitude,
        radius,
    }
}

/// Gets simplified VSOP87 data for Venus.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_venus_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Venus
    // For now, return minimal data structure
    PlanetVsop87Data {
        longitude: Vsop87Series {
            series_0: vec![],
            series_1: vec![],
            series_2: vec![],
            series_3: vec![],
            series_4: None,
            series_5: None,
        },
        latitude: Vsop87Series {
            series_0: vec![],
            series_1: vec![],
            series_2: vec![],
            series_3: vec![],
            series_4: None,
            series_5: None,
        },
        radius: Vsop87Series {
            series_0: vec![],
            series_1: vec![],
            series_2: vec![],
            series_3: vec![],
            series_4: None,
            series_5: None,
        },
    }
}

/// Gets simplified VSOP87 data for Earth.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_earth_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Earth
    get_venus_vsop87_data() // Placeholder
}

/// Gets simplified VSOP87 data for Mars.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_mars_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Mars
    get_venus_vsop87_data() // Placeholder
}

/// Gets simplified VSOP87 data for Jupiter.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_jupiter_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Jupiter
    get_venus_vsop87_data() // Placeholder
}

/// Gets simplified VSOP87 data for Saturn.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_saturn_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Saturn
    get_venus_vsop87_data() // Placeholder
}

/// Gets simplified VSOP87 data for Uranus.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_uranus_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Uranus
    get_venus_vsop87_data() // Placeholder
}

/// Gets simplified VSOP87 data for Neptune.
/// Placeholder implementation - will be populated with truncated coefficients.
fn get_neptune_vsop87_data() -> PlanetVsop87Data {
    // TODO: Add truncated VSOP87 coefficients for Neptune
    get_venus_vsop87_data() // Placeholder
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planet_from_str() {
        assert_eq!(Planet::from_str("mercury"), Some(Planet::Mercury));
        assert_eq!(Planet::from_str("JUPITER"), Some(Planet::Jupiter));
        assert_eq!(Planet::from_str("invalid"), None);
    }

    #[test]
    fn test_planet_name() {
        assert_eq!(Planet::Mars.name(), "Mars");
        assert_eq!(Planet::Saturn.name(), "Saturn");
    }

    #[test]
    fn test_evaluate_vsop87_term() {
        let term = Vsop87Term {
            amplitude: 1.0,
            phase: 0.0,
            frequency: 1.0,
        };
        // At t=0: cos(0 + 0) = 1.0
        assert!((evaluate_vsop87_term(&term, 0.0) - 1.0).abs() < 1e-10);
        
        // At t=π: cos(0 + π) = -1.0
        let t = std::f64::consts::PI;
        assert!((evaluate_vsop87_term(&term, t) + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_get_planet_vsop87_data() {
        // Test that data is available for all planets
        for planet in [Planet::Mercury, Planet::Venus, Planet::Earth, Planet::Mars,
                       Planet::Jupiter, Planet::Saturn, Planet::Uranus, Planet::Neptune] {
            let data = get_planet_vsop87_data(planet);
            assert!(data.is_some(), "VSOP87 data should be available for {}", planet.name());
        }
    }

    #[test]
    fn test_mercury_vsop87_data_structure() {
        let data = get_planet_vsop87_data(Planet::Mercury).unwrap();
        
        // Verify longitude series structure
        assert!(!data.longitude.series_0.is_empty(), "Mercury L0 series should have terms");
        assert!(!data.longitude.series_1.is_empty(), "Mercury L1 series should have terms");
        assert!(data.longitude.series_4.is_some(), "Mercury L4 series should exist");
        assert!(data.longitude.series_5.is_some(), "Mercury L5 series should exist");
        
        // Verify latitude series structure
        assert!(!data.latitude.series_0.is_empty(), "Mercury B0 series should have terms");
        assert!(data.latitude.series_4.is_none(), "Mercury B4 series should not exist");
        assert!(data.latitude.series_5.is_none(), "Mercury B5 series should not exist");
        
        // Verify radius series structure
        assert!(!data.radius.series_0.is_empty(), "Mercury R0 series should have terms");
        assert!(data.radius.series_4.is_none(), "Mercury R4 series should not exist");
        assert!(data.radius.series_5.is_none(), "Mercury R5 series should not exist");
    }

    #[test]
    fn test_vsop87_term_structure() {
        // Test that VSOP87 terms have valid structure
        let term = Vsop87Term {
            amplitude: 1.0,
            phase: 0.0,
            frequency: 26087.90314157420, // Mercury's main frequency
        };
        
        // Verify term can be evaluated
        let result = evaluate_vsop87_term(&term, 0.0);
        assert!((result - 1.0).abs() < 1e-10, "Term evaluation should work");
        
        // Verify amplitude is positive (typical for VSOP87)
        assert!(term.amplitude > 0.0, "Amplitude should typically be positive");
    }
}
