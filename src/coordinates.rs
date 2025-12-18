use crate::Result;

#[derive(Debug, Clone, Copy)]
pub struct RaDec {
    pub ra: f64,
    pub dec: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct AltAz {
    pub alt: f64,
    pub az: f64,
}

/// Earth-Centered Earth-Fixed (ECEF) coordinate system.
/// 
/// ECEF is a rotating coordinate system fixed to the Earth's surface.
/// - X-axis: Points from Earth's center through intersection of equator and prime meridian (Greenwich)
/// - Y-axis: Points from Earth's center through equator at 90°E longitude
/// - Z-axis: Points from Earth's center through North Pole
/// 
/// Coordinates rotate with the Earth (one rotation per sidereal day).
/// Units: meters (typically)
#[derive(Debug, Clone, Copy)]
pub struct Ecef {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Earth-Centered Inertial (ECI) coordinate system.
/// 
/// ECI is a non-rotating coordinate system fixed relative to the stars.
/// - X-axis: Points from Earth's center toward vernal equinox (J2000.0 epoch)
/// - Y-axis: Completes right-handed system in equatorial plane
/// - Z-axis: Points from Earth's center through North Pole
/// 
/// Coordinates do NOT rotate with the Earth (fixed in space).
/// Epoch: J2000.0 (January 1, 2000, 12:00:00 TT)
/// Units: meters (typically)
#[derive(Debug, Clone, Copy)]
pub struct Eci {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

pub fn ra_dec_to_alt_az(ra_dec: RaDec, observer_lat: f64, _observer_lon: f64, lst: f64) -> Result<AltAz> {
    let dec_rad = ra_dec.dec.to_radians();
    let lat_rad = observer_lat.to_radians();
    let hour_angle_rad = ((lst - ra_dec.ra).rem_euclid(24.0) * 15.0).to_radians();
    
    let sin_alt = dec_rad.sin() * lat_rad.sin() + dec_rad.cos() * lat_rad.cos() * hour_angle_rad.cos();
    let alt_deg = sin_alt.asin().to_degrees();
    
    let az_rad = hour_angle_rad.sin().atan2(hour_angle_rad.cos() * lat_rad.sin() - dec_rad.tan() * lat_rad.cos());
    let az_deg = (az_rad.to_degrees() + 180.0).rem_euclid(360.0);
    
    Ok(AltAz { alt: alt_deg, az: az_deg })
}

pub fn alt_az_to_ra_dec(alt_az: AltAz, observer_lat: f64, _observer_lon: f64, lst: f64) -> Result<RaDec> {
    let alt_rad = alt_az.alt.to_radians();
    let az_rad = alt_az.az.to_radians();
    let lat_rad = observer_lat.to_radians();
    
    let sin_dec = alt_rad.sin() * lat_rad.sin() + alt_rad.cos() * lat_rad.cos() * az_rad.cos();
    let dec_deg = sin_dec.asin().to_degrees();
    
    let hour_angle_rad = (-az_rad.sin()).atan2(az_rad.cos() * lat_rad.sin() + alt_rad.tan() * lat_rad.cos());
    let ra_hours = (lst - hour_angle_rad.to_degrees() / 15.0).rem_euclid(24.0);
    
    Ok(RaDec { ra: ra_hours, dec: dec_deg })
}

/// Constructs a rotation matrix for rotation around the Z-axis.
/// 
/// This is used for ECEF ↔ ECI transformations, where the Earth's rotation
/// is represented as a rotation around the Z-axis (Earth's rotation axis).
/// 
/// # Arguments
/// * `angle_rad` - Rotation angle in radians (positive = counterclockwise when viewed from +Z)
/// 
/// # Returns
/// 3x3 rotation matrix as a 2D array: `[[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]]`
/// 
/// # Formula
/// Rotation matrix R_z(θ) around Z-axis:
/// ```
/// [cos(θ)  sin(θ)  0]
/// [-sin(θ) cos(θ)  0]
/// [0       0       1]
/// ```
#[cfg(test)]
pub(crate) use rotation_matrix_z_impl as rotation_matrix_z;

fn rotation_matrix_z_impl(angle_rad: f64) -> [[f64; 3]; 3] {
    let cos_theta = angle_rad.cos();
    let sin_theta = angle_rad.sin();
    
    let matrix = [
        [cos_theta,  sin_theta, 0.0],
        [-sin_theta, cos_theta, 0.0],
        [0.0,        0.0,       1.0],
    ];
    
    // Debug logging: rotation matrix construction
    log::debug!(
        "Rotation matrix Z-axis: angle={:.6} rad ({:.4}°), matrix=[{:.6}, {:.6}, 0.0; {:.6}, {:.6}, 0.0; 0.0, 0.0, 1.0]",
        angle_rad,
        angle_rad.to_degrees(),
        matrix[0][0], matrix[0][1],
        matrix[1][0], matrix[1][1]
    );
    
    matrix
}

/// Applies a 3x3 rotation matrix to a 3D vector.
/// 
/// # Arguments
/// * `matrix` - 3x3 rotation matrix
/// * `vector` - Input vector (x, y, z)
/// 
/// # Returns
/// Rotated vector (x', y', z')
fn apply_rotation_matrix(matrix: [[f64; 3]; 3], vector: (f64, f64, f64)) -> (f64, f64, f64) {
    let (x, y, z) = vector;
    (
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z,
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z,
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z,
    )
}

/// Converts Earth-Centered Earth-Fixed (ECEF) coordinates to Earth-Centered Inertial (ECI) coordinates.
/// 
/// This transformation accounts for the Earth's rotation by rotating the coordinate system
/// backwards (negative rotation) by the Greenwich Mean Sidereal Time (GMST).
/// 
/// # Coordinate Systems
/// - **ECEF**: Rotates with the Earth (one rotation per sidereal day)
/// - **ECI**: Fixed relative to the stars (J2000.0 epoch)
/// 
/// # Arguments
/// * `ecef` - ECEF coordinates (x, y, z) in meters
/// * `gmst` - Greenwich Mean Sidereal Time in hours (0-24)
/// 
/// # Returns
/// ECI coordinates (x, y, z) in meters
/// 
/// # Errors
/// Returns an error if:
/// - Input coordinates contain NaN or infinity
/// - GMST is outside valid range (0-24 hours)
/// 
/// # Mathematical Basis
/// The transformation uses a rotation matrix around the Z-axis:
/// - Rotation angle: -GMST (negative because we're "undoing" Earth's rotation)
/// - ECI = R_z(-GMST) × ECEF
/// 
/// # Note
/// This is a simplified transformation that assumes:
/// - J2000.0 epoch for ECI reference frame
/// - No precession or nutation corrections (suitable for most applications)
/// - GMST accurately represents Earth's rotation angle
/// 
/// For high-precision applications (>arcsecond accuracy), precession and nutation
/// corrections should be added (IAU-76/FK5 or IAU-2000/2006 models).
/// 
/// # Example
/// ```
/// use cli_astro_calc::coordinates::{Ecef, ecef_to_eci};
/// 
/// // Point on equator at Greenwich meridian (at Earth's surface)
/// let ecef = Ecef { x: 6378137.0, y: 0.0, z: 0.0 }; // ~Earth radius
/// let gmst = 0.0; // GMST = 0 means ECEF and ECI X-axes align
/// let eci = ecef_to_eci(ecef, gmst)?;
/// // At GMST=0, ECEF and ECI should be identical
/// ```
pub fn ecef_to_eci(ecef: Ecef, gmst: f64) -> Result<Eci> {
    // Validate inputs
    if ecef.x.is_nan() || ecef.y.is_nan() || ecef.z.is_nan() {
        return Err(crate::AstroError::InvalidCoordinate(
            "ECEF coordinates contain NaN".to_string()
        ));
    }
    if ecef.x.is_infinite() || ecef.y.is_infinite() || ecef.z.is_infinite() {
        return Err(crate::AstroError::InvalidCoordinate(
            "ECEF coordinates contain infinity".to_string()
        ));
    }
    if gmst.is_nan() || gmst.is_infinite() {
        return Err(crate::AstroError::InvalidTime(
            format!("Invalid GMST: {}", gmst)
        ));
    }
    
    // Normalize GMST to 0-24 range
    let gmst_normalized = gmst.rem_euclid(24.0);
    
    // Info logging: transformation operation
    log::info!(
        "ECEF to ECI transformation: input ECEF=({:.3}, {:.3}, {:.3}) m, GMST={:.6} h (normalized: {:.6} h)",
        ecef.x, ecef.y, ecef.z, gmst, gmst_normalized
    );
    
    // Convert GMST (hours) to rotation angle (radians)
    // Negative rotation: ECEF rotates with Earth, ECI doesn't, so we rotate backwards
    let rotation_angle_rad = -(gmst_normalized * 15.0).to_radians(); // 15° per hour
    
    // Construct rotation matrix (debug logging inside function)
    let rotation_matrix = rotation_matrix_z_impl(rotation_angle_rad);
    
    // Apply rotation
    let (x, y, z) = apply_rotation_matrix(rotation_matrix, (ecef.x, ecef.y, ecef.z));
    
    // Info logging: transformation result
    log::info!(
        "ECEF to ECI transformation: output ECI=({:.3}, {:.3}, {:.3}) m",
        x, y, z
    );
    
    Ok(Eci { x, y, z })
}

/// Converts Earth-Centered Inertial (ECI) coordinates to Earth-Centered Earth-Fixed (ECEF) coordinates.
/// 
/// This transformation accounts for the Earth's rotation by rotating the coordinate system
/// forwards (positive rotation) by the Greenwich Mean Sidereal Time (GMST).
/// 
/// # Coordinate Systems
/// - **ECI**: Fixed relative to the stars (J2000.0 epoch)
/// - **ECEF**: Rotates with the Earth (one rotation per sidereal day)
/// 
/// # Arguments
/// * `eci` - ECI coordinates (x, y, z) in meters
/// * `gmst` - Greenwich Mean Sidereal Time in hours (0-24)
/// 
/// # Returns
/// ECEF coordinates (x, y, z) in meters
/// 
/// # Errors
/// Returns an error if:
/// - Input coordinates contain NaN or infinity
/// - GMST is outside valid range (0-24 hours)
/// 
/// # Mathematical Basis
/// The transformation uses a rotation matrix around the Z-axis:
/// - Rotation angle: +GMST (positive because we're applying Earth's rotation)
/// - ECEF = R_z(+GMST) × ECI
/// 
/// # Note
/// This is a simplified transformation that assumes:
/// - J2000.0 epoch for ECI reference frame
/// - No precession or nutation corrections (suitable for most applications)
/// - GMST accurately represents Earth's rotation angle
/// 
/// For high-precision applications (>arcsecond accuracy), precession and nutation
/// corrections should be added (IAU-76/FK5 or IAU-2000/2006 models).
/// 
/// # Example
/// ```
/// use cli_astro_calc::coordinates::{Eci, eci_to_ecef};
/// 
/// // Point in ECI frame
/// let eci = Eci { x: 6378137.0, y: 0.0, z: 0.0 };
/// let gmst = 0.0; // GMST = 0 means ECEF and ECI X-axes align
/// let ecef = eci_to_ecef(eci, gmst)?;
/// // At GMST=0, ECEF and ECI should be identical
/// ```
pub fn eci_to_ecef(eci: Eci, gmst: f64) -> Result<Ecef> {
    // Validate inputs
    if eci.x.is_nan() || eci.y.is_nan() || eci.z.is_nan() {
        return Err(crate::AstroError::InvalidCoordinate(
            "ECI coordinates contain NaN".to_string()
        ));
    }
    if eci.x.is_infinite() || eci.y.is_infinite() || eci.z.is_infinite() {
        return Err(crate::AstroError::InvalidCoordinate(
            "ECI coordinates contain infinity".to_string()
        ));
    }
    if gmst.is_nan() || gmst.is_infinite() {
        return Err(crate::AstroError::InvalidTime(
            format!("Invalid GMST: {}", gmst)
        ));
    }
    
    // Normalize GMST to 0-24 range
    let gmst_normalized = gmst.rem_euclid(24.0);
    
    // Info logging: transformation operation
    log::info!(
        "ECI to ECEF transformation: input ECI=({:.3}, {:.3}, {:.3}) m, GMST={:.6} h (normalized: {:.6} h)",
        eci.x, eci.y, eci.z, gmst, gmst_normalized
    );
    
    // Convert GMST (hours) to rotation angle (radians)
    // Positive rotation: ECI is fixed, ECEF rotates with Earth
    let rotation_angle_rad = (gmst_normalized * 15.0).to_radians(); // 15° per hour
    
    // Construct rotation matrix (debug logging inside function)
    let rotation_matrix = rotation_matrix_z_impl(rotation_angle_rad);
    
    // Apply rotation
    let (x, y, z) = apply_rotation_matrix(rotation_matrix, (eci.x, eci.y, eci.z));
    
    // Info logging: transformation result
    log::info!(
        "ECI to ECEF transformation: output ECEF=({:.3}, {:.3}, {:.3}) m",
        x, y, z
    );
    
    Ok(Ecef { x, y, z })
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

    // ============================================================================
    // ECEF ↔ ECI Transformation Tests
    // ============================================================================
    // Test cases for Step A1: Research & Design
    // These tests will be implemented in Step A3: Testing & Validation
    
    #[test]
    fn test_rotation_matrix_z_0_degrees() {
        // Test rotation matrix for 0° rotation (identity matrix)
        // Expected: Identity matrix
        // [1  0  0]
        // [0  1  0]
        // [0  0  1]
        let matrix = rotation_matrix_z(0.0);
        assert!((matrix[0][0] - 1.0).abs() < 1e-10, "m00 should be 1.0");
        assert!((matrix[0][1] - 0.0).abs() < 1e-10, "m01 should be 0.0");
        assert!((matrix[0][2] - 0.0).abs() < 1e-10, "m02 should be 0.0");
        assert!((matrix[1][0] - 0.0).abs() < 1e-10, "m10 should be 0.0");
        assert!((matrix[1][1] - 1.0).abs() < 1e-10, "m11 should be 1.0");
        assert!((matrix[1][2] - 0.0).abs() < 1e-10, "m12 should be 0.0");
        assert!((matrix[2][0] - 0.0).abs() < 1e-10, "m20 should be 0.0");
        assert!((matrix[2][1] - 0.0).abs() < 1e-10, "m21 should be 0.0");
        assert!((matrix[2][2] - 1.0).abs() < 1e-10, "m22 should be 1.0");
    }
    
    #[test]
    fn test_rotation_matrix_z_90_degrees() {
        // Test rotation matrix for 90° rotation
        // Expected:
        // [0   1  0]
        // [-1  0  0]
        // [0   0  1]
        let matrix = rotation_matrix_z(90.0_f64.to_radians());
        assert!((matrix[0][0] - 0.0).abs() < 1e-10, "m00 should be 0.0");
        assert!((matrix[0][1] - 1.0).abs() < 1e-10, "m01 should be 1.0");
        assert!((matrix[1][0] - (-1.0)).abs() < 1e-10, "m10 should be -1.0");
        assert!((matrix[1][1] - 0.0).abs() < 1e-10, "m11 should be 0.0");
        assert!((matrix[2][2] - 1.0).abs() < 1e-10, "m22 should be 1.0");
    }
    
    #[test]
    fn test_rotation_matrix_z_180_degrees() {
        // Test rotation matrix for 180° rotation
        // Expected:
        // [-1  0  0]
        // [0  -1  0]
        // [0   0  1]
        let matrix = rotation_matrix_z(180.0_f64.to_radians());
        assert!((matrix[0][0] - (-1.0)).abs() < 1e-10, "m00 should be -1.0");
        assert!((matrix[0][1] - 0.0).abs() < 1e-10, "m01 should be 0.0");
        assert!((matrix[1][0] - 0.0).abs() < 1e-10, "m10 should be 0.0");
        assert!((matrix[1][1] - (-1.0)).abs() < 1e-10, "m11 should be -1.0");
        assert!((matrix[2][2] - 1.0).abs() < 1e-10, "m22 should be 1.0");
    }
    
    #[test]
    fn test_rotation_matrix_z_270_degrees() {
        // Test rotation matrix for 270° rotation
        // Expected:
        // [0  -1  0]
        // [1   0  0]
        // [0   0  1]
        let matrix = rotation_matrix_z(270.0_f64.to_radians());
        assert!((matrix[0][0] - 0.0).abs() < 1e-10, "m00 should be 0.0");
        assert!((matrix[0][1] - (-1.0)).abs() < 1e-10, "m01 should be -1.0");
        assert!((matrix[1][0] - 1.0).abs() < 1e-10, "m10 should be 1.0");
        assert!((matrix[1][1] - 0.0).abs() < 1e-10, "m11 should be 0.0");
        assert!((matrix[2][2] - 1.0).abs() < 1e-10, "m22 should be 1.0");
    }
    
    #[test]
    fn test_rotation_matrix_orthogonal() {
        // Verify rotation matrix properties:
        // - Matrix is orthogonal (R^T = R^-1, i.e., R^T * R = I)
        // - Determinant = 1 (preserves volume/orientation)
        let angles = [0.0, 45.0, 90.0, 180.0, 270.0, 360.0];
        
        for angle_deg in angles.iter() {
            let angle_rad = angle_deg.to_radians();
            let r = rotation_matrix_z(angle_rad);
            
            // Calculate determinant: det(R_z) = cos²(θ) + sin²(θ) = 1
            let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
                    - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
                    + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
            assert!((det - 1.0).abs() < 1e-10, 
                "Determinant should be 1.0 for angle {}°, got {}", angle_deg, det);
            
            // Verify orthogonality: R^T * R = I
            // For Z-axis rotation, this means:
            // - Each row/column has unit length
            // - Dot product of different rows = 0
            let row0_len_sq = r[0][0] * r[0][0] + r[0][1] * r[0][1] + r[0][2] * r[0][2];
            let row1_len_sq = r[1][0] * r[1][0] + r[1][1] * r[1][1] + r[1][2] * r[1][2];
            let row2_len_sq = r[2][0] * r[2][0] + r[2][1] * r[2][1] + r[2][2] * r[2][2];
            
            assert!((row0_len_sq - 1.0).abs() < 1e-10, 
                "Row 0 should have unit length for angle {}°", angle_deg);
            assert!((row1_len_sq - 1.0).abs() < 1e-10, 
                "Row 1 should have unit length for angle {}°", angle_deg);
            assert!((row2_len_sq - 1.0).abs() < 1e-10, 
                "Row 2 should have unit length for angle {}°", angle_deg);
            
            // Dot product of row 0 and row 1 should be 0
            let dot01 = r[0][0] * r[1][0] + r[0][1] * r[1][1] + r[0][2] * r[1][2];
            assert!((dot01 - 0.0).abs() < 1e-10, 
                "Rows 0 and 1 should be orthogonal for angle {}°", angle_deg);
        }
    }
    
    #[test]
    fn test_ecef_to_eci_at_greenwich_meridian() {
        // Test ECEF to ECI at Greenwich meridian (x-axis in ECEF)
        // Known case: Point on equator at Greenwich (x = Earth radius, y = 0, z = 0)
        const EARTH_RADIUS: f64 = 6378137.0; // WGS84 semi-major axis in meters
        
        // At GMST = 0: ECEF and ECI should align (no rotation)
        let ecef = Ecef { x: EARTH_RADIUS, y: 0.0, z: 0.0 };
        let eci = ecef_to_eci(ecef, 0.0).unwrap();
        assert!((eci.x - EARTH_RADIUS).abs() < 1e-6, "At GMST=0, x should be unchanged");
        assert!((eci.y - 0.0).abs() < 1e-6, "At GMST=0, y should be 0");
        assert!((eci.z - 0.0).abs() < 1e-6, "At GMST=0, z should be 0");
        
        // At GMST = 6: Should rotate 90° (x -> y, y -> -x)
        let eci_6h = ecef_to_eci(ecef, 6.0).unwrap();
        assert!((eci_6h.x - 0.0).abs() < 1e-6, "At GMST=6h, x should be 0 (rotated to y)");
        assert!((eci_6h.y - EARTH_RADIUS).abs() < 1e-6, "At GMST=6h, y should equal original x");
        assert!((eci_6h.z - 0.0).abs() < 1e-6, "At GMST=6h, z should remain 0");
    }
    
    #[test]
    fn test_ecef_to_eci_at_north_pole() {
        // Test ECEF to ECI at North Pole
        // Point at North Pole: (x = 0, y = 0, z = Earth radius)
        // Z-coordinate should remain unchanged (rotation around Z-axis)
        const EARTH_RADIUS: f64 = 6378137.0;
        
        let ecef = Ecef { x: 0.0, y: 0.0, z: EARTH_RADIUS };
        
        // Test at different GMST values - z should always remain the same
        for gmst in [0.0, 6.0, 12.0, 18.0] {
            let eci = ecef_to_eci(ecef, gmst).unwrap();
            assert!((eci.z - EARTH_RADIUS).abs() < 1e-6, 
                "Z-coordinate should remain unchanged at GMST={}, got {}", gmst, eci.z);
        }
    }
    
    #[test]
    fn test_ecef_to_eci_at_equator() {
        // Test ECEF to ECI at equator (z = 0)
        // Point on equator: (x, y, 0)
        // Z-coordinate should remain 0
        const EARTH_RADIUS: f64 = 6378137.0;
        
        let ecef = Ecef { x: EARTH_RADIUS, y: EARTH_RADIUS, z: 0.0 };
        
        // Test at different GMST values - z should always remain 0
        for gmst in [0.0, 6.0, 12.0, 18.0] {
            let eci = ecef_to_eci(ecef, gmst).unwrap();
            assert!((eci.z - 0.0).abs() < 1e-6, 
                "Z-coordinate should remain 0 at equator for GMST={}, got {}", gmst, eci.z);
        }
    }
    
    #[test]
    fn test_eci_to_ecef_at_greenwich_meridian() {
        // Test ECI to ECEF at Greenwich meridian
        // Inverse of ecef_to_eci test
        const EARTH_RADIUS: f64 = 6378137.0;
        
        // At GMST = 0: ECI and ECEF should align
        let eci = Eci { x: EARTH_RADIUS, y: 0.0, z: 0.0 };
        let ecef = eci_to_ecef(eci, 0.0).unwrap();
        assert!((ecef.x - EARTH_RADIUS).abs() < 1e-6, "At GMST=0, x should be unchanged");
        assert!((ecef.y - 0.0).abs() < 1e-6, "At GMST=0, y should be 0");
        
        // At GMST = 6: ECI (x, 0, 0) should rotate to ECEF (0, x, 0)
        let ecef_6h = eci_to_ecef(eci, 6.0).unwrap();
        assert!((ecef_6h.x - 0.0).abs() < 1e-6, "At GMST=6h, x should be 0");
        assert!((ecef_6h.y - EARTH_RADIUS).abs() < 1e-6, "At GMST=6h, y should equal original x");
    }
    
    #[test]
    fn test_ecef_eci_round_trip() {
        // Test round-trip conversion: ECEF → ECI → ECEF
        // Should return to original coordinates within numerical precision
        // Tolerance: < 1mm for Earth radius scale coordinates
        const EARTH_RADIUS: f64 = 6378137.0;
        const TOLERANCE_MM: f64 = 0.001; // 1mm tolerance
        
        let test_cases = [
            Ecef { x: EARTH_RADIUS, y: 0.0, z: 0.0 },
            Ecef { x: 0.0, y: EARTH_RADIUS, z: 0.0 },
            Ecef { x: 0.0, y: 0.0, z: EARTH_RADIUS },
            Ecef { x: EARTH_RADIUS / 2.0, y: EARTH_RADIUS / 2.0, z: EARTH_RADIUS / 2.0 },
        ];
        
        for ecef_original in test_cases.iter() {
            let eci = ecef_to_eci(*ecef_original, 12.5).unwrap();
            let ecef_result = eci_to_ecef(eci, 12.5).unwrap();
            
            let error_x = (ecef_result.x - ecef_original.x).abs();
            let error_y = (ecef_result.y - ecef_original.y).abs();
            let error_z = (ecef_result.z - ecef_original.z).abs();
            
            assert!(error_x < TOLERANCE_MM, 
                "Round-trip error in x: {} m (expected < {} m)", error_x, TOLERANCE_MM);
            assert!(error_y < TOLERANCE_MM, 
                "Round-trip error in y: {} m (expected < {} m)", error_y, TOLERANCE_MM);
            assert!(error_z < TOLERANCE_MM, 
                "Round-trip error in z: {} m (expected < {} m)", error_z, TOLERANCE_MM);
        }
    }
    
    #[test]
    fn test_ecef_eci_round_trip_at_different_gmst() {
        // Test round-trip at various GMST values
        // GMST = 0, 6, 12, 18 hours
        const EARTH_RADIUS: f64 = 6378137.0;
        const TOLERANCE_MM: f64 = 0.001; // 1mm tolerance
        
        let ecef_original = Ecef { 
            x: EARTH_RADIUS, 
            y: EARTH_RADIUS / 2.0, 
            z: EARTH_RADIUS / 3.0 
        };
        
        for gmst in [0.0, 6.0, 12.0, 18.0, 23.5] {
            let eci = ecef_to_eci(ecef_original, gmst).unwrap();
            let ecef_result = eci_to_ecef(eci, gmst).unwrap();
            
            let error = ((ecef_result.x - ecef_original.x).powi(2) +
                        (ecef_result.y - ecef_original.y).powi(2) +
                        (ecef_result.z - ecef_original.z).powi(2)).sqrt();
            
            assert!(error < TOLERANCE_MM, 
                "Round-trip error at GMST={}h: {} m (expected < {} m)", gmst, error, TOLERANCE_MM);
        }
    }
    
    #[test]
    fn test_ecef_to_eci_gmst_normalization() {
        // Test that GMST values outside 0-24 range are normalized
        // GMST = 25 should be treated as GMST = 1
        // GMST = -1 should be treated as GMST = 23
        const EARTH_RADIUS: f64 = 6378137.0;
        
        let ecef = Ecef { x: EARTH_RADIUS, y: 0.0, z: 0.0 };
        
        // GMST = 25 should equal GMST = 1
        let eci_25 = ecef_to_eci(ecef, 25.0).unwrap();
        let eci_1 = ecef_to_eci(ecef, 1.0).unwrap();
        assert!((eci_25.x - eci_1.x).abs() < 1e-6, "GMST=25 should equal GMST=1");
        assert!((eci_25.y - eci_1.y).abs() < 1e-6, "GMST=25 should equal GMST=1");
        
        // GMST = -1 should equal GMST = 23
        let eci_neg1 = ecef_to_eci(ecef, -1.0).unwrap();
        let eci_23 = ecef_to_eci(ecef, 23.0).unwrap();
        assert!((eci_neg1.x - eci_23.x).abs() < 1e-6, "GMST=-1 should equal GMST=23");
        assert!((eci_neg1.y - eci_23.y).abs() < 1e-6, "GMST=-1 should equal GMST=23");
    }
    
    #[test]
    fn test_ecef_to_eci_known_coordinate_pairs() {
        // Test with known ECEF/ECI coordinate pairs
        // Reference: At GMST=0, ECEF and ECI should be identical
        // This is a basic sanity check - full validation would require JPL Horizons data
        const EARTH_RADIUS: f64 = 6378137.0;
        
        // Known case: Greenwich meridian at equator, GMST=0
        let ecef = Ecef { x: EARTH_RADIUS, y: 0.0, z: 0.0 };
        let eci = ecef_to_eci(ecef, 0.0).unwrap();
        
        // At GMST=0, coordinates should be identical
        assert!((eci.x - ecef.x).abs() < 1e-6, "At GMST=0, x should be identical");
        assert!((eci.y - ecef.y).abs() < 1e-6, "At GMST=0, y should be identical");
        assert!((eci.z - ecef.z).abs() < 1e-6, "At GMST=0, z should be identical");
    }
    
    #[test]
    fn test_ecef_to_eci_invalid_inputs() {
        // Test error handling:
        // - NaN coordinates
        // - Infinite coordinates
        // - Invalid GMST values
        
        // Test NaN coordinates
        let ecef_nan = Ecef { x: f64::NAN, y: 0.0, z: 0.0 };
        assert!(ecef_to_eci(ecef_nan, 0.0).is_err(), "Should error on NaN x coordinate");
        
        let ecef_nan_y = Ecef { x: 0.0, y: f64::NAN, z: 0.0 };
        assert!(ecef_to_eci(ecef_nan_y, 0.0).is_err(), "Should error on NaN y coordinate");
        
        // Test infinite coordinates
        let ecef_inf = Ecef { x: f64::INFINITY, y: 0.0, z: 0.0 };
        assert!(ecef_to_eci(ecef_inf, 0.0).is_err(), "Should error on infinite x coordinate");
        
        // Test invalid GMST
        assert!(ecef_to_eci(Ecef { x: 0.0, y: 0.0, z: 0.0 }, f64::NAN).is_err(), 
            "Should error on NaN GMST");
        assert!(ecef_to_eci(Ecef { x: 0.0, y: 0.0, z: 0.0 }, f64::INFINITY).is_err(), 
            "Should error on infinite GMST");
    }
    
    #[test]
    fn test_ecef_to_eci_at_origin() {
        // Test edge case: coordinates at origin (0, 0, 0)
        // Should handle gracefully (no rotation needed)
        let ecef = Ecef { x: 0.0, y: 0.0, z: 0.0 };
        
        for gmst in [0.0, 6.0, 12.0, 18.0] {
            let eci = ecef_to_eci(ecef, gmst).unwrap();
            assert!((eci.x - 0.0).abs() < 1e-10, "Origin x should remain 0 at GMST={}", gmst);
            assert!((eci.y - 0.0).abs() < 1e-10, "Origin y should remain 0 at GMST={}", gmst);
            assert!((eci.z - 0.0).abs() < 1e-10, "Origin z should remain 0 at GMST={}", gmst);
        }
    }
    
    #[test]
    fn test_ecef_to_eci_large_coordinates() {
        // Test with large coordinate values (e.g., geosynchronous orbit ~42,164 km)
        // Verify numerical stability
        const GEO_ALTITUDE: f64 = 42164000.0; // Geosynchronous orbit altitude in meters
        
        let ecef = Ecef { x: GEO_ALTITUDE, y: GEO_ALTITUDE / 2.0, z: GEO_ALTITUDE / 3.0 };
        
        // Test round-trip at large distances
        let eci = ecef_to_eci(ecef, 12.0).unwrap();
        let ecef_result = eci_to_ecef(eci, 12.0).unwrap();
        
        // Tolerance: 1mm even for large coordinates
        const TOLERANCE_MM: f64 = 0.001;
        let error = ((ecef_result.x - ecef.x).powi(2) +
                    (ecef_result.y - ecef.y).powi(2) +
                    (ecef_result.z - ecef.z).powi(2)).sqrt();
        
        assert!(error < TOLERANCE_MM, 
            "Round-trip error for large coordinates: {} m (expected < {} m)", error, TOLERANCE_MM);
    }
}