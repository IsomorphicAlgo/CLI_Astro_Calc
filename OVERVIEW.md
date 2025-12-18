# CLI Astro Calc - Overview

## Core Purpose

**CLI Astro Calc** is a command-line astronomy and orbital mechanics calculator built in Rust. The program provides accurate astronomical calculations for both educational and practical applications, including:

- **Time Systems**: Converting between calendar dates and astronomical time systems (Julian Date, Sidereal Time)
- **Celestial Object Positions**: Calculating the positions of the Sun and Moon in the sky
- **Rise/Set Times**: Determining when celestial objects rise and set for any location on Earth
- **Coordinate Transformations**: Converting between different astronomical coordinate systems
- **Orbital Mechanics**: Solving Kepler's equation and converting between orbital elements and state vectors
- **Earth-Centered Coordinates**: Transforming between Earth-Centered Earth-Fixed (ECEF) and Earth-Centered Inertial (ECI) coordinate systems

The program is designed for accuracy, with all calculations verified against authoritative astronomical sources and comprehensive test coverage.

---

## Top-Level Feature Summary

### 1. Time Systems (`time.rs`)
- **Julian Date Calculation**: Converts calendar dates to continuous day count since 4713 BC
- **Greenwich Mean Sidereal Time (GMST)**: Calculates Earth's rotation relative to the stars
- **Local Sidereal Time (LST)**: Converts GMST to local sidereal time for any longitude

### 2. Celestial Object Positions (`celestial.rs`)
- **Solar Position**: Calculates Sun's Right Ascension (RA) and Declination (Dec) using mean anomaly and equation of center
- **Lunar Position**: Calculates Moon's RA/Dec using perturbation theory with multiple periodic terms
- **Rise/Set Times**: Determines when Sun or Moon rises/sets for any observer location

### 3. Coordinate Conversions (`coordinates.rs`)
- **RA/Dec ↔ Alt/Az**: Converts between equatorial (celestial) and horizontal (observer-based) coordinate systems
- **ECEF ↔ ECI**: Transforms between Earth-Centered Earth-Fixed and Earth-Centered Inertial coordinate systems using rotation matrices

### 4. Orbital Mechanics (`orbital.rs`)
- **Kepler's Equation Solver**: Iteratively solves for eccentric anomaly from mean anomaly
- **Orbital Period**: Calculates orbital period using Kepler's Third Law
- **Elements to State Vectors**: Converts classical orbital elements to position/velocity vectors

---

## Mathematical Equations and Their Applications

### Time Systems

#### Julian Date (JD)
**Purpose**: Provides a continuous day count for astronomical calculations, independent of calendar systems.

**Formula**:
```
JD = day + (153×m + 2)/5 + 365×y + y/4 - y/100 + y/400 - 32045 + time_fraction
```
Where:
- `m` = month (adjusted for algorithm)
- `y` = year (adjusted for algorithm)
- `time_fraction` = (hour - 12 + minute/60 + second/3600) / 24

**Use**: Foundation for all time-dependent astronomical calculations. J2000.0 epoch = JD 2451545.0.

#### Greenwich Mean Sidereal Time (GMST)
**Purpose**: Measures Earth's rotation angle relative to the vernal equinox (fixed stars), not the Sun.

**Formula**:
```
GMST = (GMST_J2000 + SIDEREAL_RATE × (JD - J2000)) mod 24
```
Where:
- `GMST_J2000` = 18.697374558 hours (GMST at J2000.0)
- `SIDEREAL_RATE` = 24.06570982441908 hours/day (sidereal day rate)

**Use**: Required for coordinate transformations between rotating (ECEF) and non-rotating (ECI) frames, and for converting between equatorial and horizontal coordinates.

#### Local Sidereal Time (LST)
**Purpose**: Converts GMST to local sidereal time for a specific observer longitude.

**Formula**:
```
LST = (GMST + longitude/15) mod 24
```
Where longitude is in degrees (positive = East, negative = West).

**Use**: Essential for converting between RA/Dec and Alt/Az coordinates, as it relates the observer's local meridian to the celestial coordinate system.

---

### Celestial Object Positions

#### Solar Position (Simplified Model)
**Purpose**: Calculates the Sun's position in equatorial coordinates (RA/Dec) for any date.

**Key Equations**:
1. **Mean Anomaly**:
   ```
   M = 357.5291° + 0.98560028° × d
   ```
   Where `d` = days since J2000.0

2. **Equation of Center** (correction for elliptical orbit):
   ```
   EOC = 1.9148° × sin(M) + 0.0200° × sin(2M) + 0.0003° × sin(3M)
   ```

3. **Ecliptic Longitude**:
   ```
   λ = M + EOC + 180° + 102.9372°
   ```

4. **Obliquity of Ecliptic**:
   ```
   ε = 23.4393° - 0.0000004° × d
   ```

5. **Equatorial Coordinates** (ecliptic → equatorial):
   ```
   RA = atan2(sin(λ) × cos(ε), cos(λ))
   Dec = arcsin(sin(λ) × sin(ε))
   ```

**Use**: Determines where the Sun appears in the sky, essential for sunrise/sunset calculations and solar observations.

#### Lunar Position (Perturbation Theory)
**Purpose**: Calculates the Moon's position using multiple periodic terms to account for its complex orbit.

**Key Equations**:
1. **Mean Arguments** (multiple angles):
   - Mean longitude: `L' = 218.316° + 481267.881° × t + ...`
   - Mean elongation: `D = 297.850° + 445267.111° × t + ...`
   - Mean anomaly (Sun): `M = 357.529° + 35999.050° × t + ...`
   - Mean anomaly (Moon): `M' = 134.963° + 477198.868° × t + ...`
   - Argument of latitude: `F = 93.272° + 483202.018° × t + ...`

2. **Perturbation Series** (sum of periodic terms):
   ```
   σ_L = 6288774" × sin(M') + 1274027" × sin(2D - M') + ...
   σ_B = 5128122" × sin(F) + 280602" × sin(M' + F) + ...
   ```

3. **Ecliptic Coordinates**:
   ```
   λ = (L' + σ_L/1000000) mod 360°
   β = (σ_B/1000000) mod 360°
   ```

4. **Equatorial Coordinates** (conversion from ecliptic):
   ```
   RA = atan2(sin(λ) × cos(ε) - tan(β) × sin(ε), cos(λ))
   Dec = arcsin(sin(β) × cos(ε) + cos(β) × sin(ε) × sin(λ))
   ```

**Use**: Accurate lunar position calculation accounting for perturbations from the Sun and Earth's non-spherical shape.

---

### Coordinate Transformations

#### RA/Dec to Alt/Az (Equatorial to Horizontal)
**Purpose**: Converts celestial coordinates (fixed relative to stars) to observer-based coordinates (altitude/azimuth).

**Key Equations**:
1. **Hour Angle**:
   ```
   HA = LST - RA (in hours, converted to degrees: HA° = HA × 15)
   ```

2. **Altitude**:
   ```
   sin(Alt) = sin(Dec) × sin(Lat) + cos(Dec) × cos(Lat) × cos(HA)
   Alt = arcsin(sin(Alt))
   ```

3. **Azimuth**:
   ```
   Az = atan2(sin(HA), cos(HA) × sin(Lat) - tan(Dec) × cos(Lat)) + 180°
   Az = Az mod 360°
   ```

**Use**: Determines where a celestial object appears in the sky from a specific location on Earth (e.g., "30° above the horizon, 180° azimuth = due South").

#### Alt/Az to RA/Dec (Horizontal to Equatorial)
**Purpose**: Inverse transformation - converts observer-based coordinates to celestial coordinates.

**Key Equations**:
1. **Declination**:
   ```
   sin(Dec) = sin(Alt) × sin(Lat) + cos(Alt) × cos(Lat) × cos(Az)
   Dec = arcsin(sin(Dec))
   ```

2. **Hour Angle**:
   ```
   HA = atan2(-sin(Az), cos(Az) × sin(Lat) + tan(Alt) × cos(Lat))
   ```

3. **Right Ascension**:
   ```
   RA = (LST - HA/15) mod 24 hours
   ```

**Use**: Converts telescope pointing coordinates (Alt/Az) to celestial coordinates for catalog lookup or tracking.

---

### ECEF/ECI Transformations

#### Rotation Matrix (Z-axis)
**Purpose**: Rotates coordinates around the Z-axis (Earth's rotation axis) to account for Earth's rotation.

**Formula**:
```
R_z(θ) = [cos(θ)   sin(θ)   0 ]
         [-sin(θ)  cos(θ)    0 ]
         [0        0         1 ]
```

**Use**: Core transformation for converting between rotating (ECEF) and non-rotating (ECI) coordinate systems.

#### ECEF to ECI
**Purpose**: Converts from Earth-Centered Earth-Fixed (rotates with Earth) to Earth-Centered Inertial (fixed relative to stars).

**Key Equations**:
1. **Rotation Angle**:
   ```
   θ = -GMST × 15° (negative because we "undo" Earth's rotation)
   ```

2. **Matrix Application**:
   ```
   [ECI_x]   [cos(θ)   sin(θ)   0 ] [ECEF_x]
   [ECI_y] = [-sin(θ)  cos(θ)    0 ] [ECEF_y]
   [ECI_z]   [0        0         1 ] [ECEF_z]
   ```

**Use**: Essential for satellite tracking, GPS calculations, and any application requiring coordinates fixed relative to stars rather than Earth's surface.

#### ECI to ECEF
**Purpose**: Inverse transformation - converts from inertial to Earth-fixed coordinates.

**Key Equations**:
1. **Rotation Angle**:
   ```
   θ = +GMST × 15° (positive because we apply Earth's rotation)
   ```

2. **Matrix Application**: Same rotation matrix as above, but with positive angle.

**Use**: Converting satellite positions from inertial frame to Earth-fixed frame for ground station tracking.

---

### Orbital Mechanics

#### Kepler's Equation
**Purpose**: Relates mean anomaly (time-based) to eccentric anomaly (geometry-based) for elliptical orbits.

**Equation**:
```
M = E - e × sin(E)
```
Where:
- `M` = Mean anomaly (degrees)
- `E` = Eccentric anomaly (degrees)
- `e` = Orbital eccentricity (0-1)

**Solution Method**: Iterative Newton-Raphson method:
```
E_{n+1} = E_n - (E_n - e × sin(E_n) - M) / (1 - e × cos(E_n))
```

**Use**: Required to determine an object's position in its orbit at a given time.

#### True Anomaly from Eccentric Anomaly
**Purpose**: Converts eccentric anomaly to true anomaly (actual angle from periapsis).

**Formula**:
```
ν = 2 × arctan(√((1+e)/(1-e)) × tan(E/2))
```

**Use**: Determines the actual angular position of an object in its orbit.

#### Orbital Period (Kepler's Third Law)
**Purpose**: Calculates the time for one complete orbit.

**Formula**:
```
T = 2π × √(a³ / μ)
```
Where:
- `T` = Orbital period (seconds)
- `a` = Semi-major axis (meters)
- `μ` = Standard gravitational parameter (m³/s²)

**Use**: Determines orbital period for satellites, planets, or any orbiting body.

#### Orbital Elements to State Vectors
**Purpose**: Converts classical orbital elements (a, e, i, Ω, ω, M) to position and velocity vectors.

**Key Steps**:
1. Solve Kepler's equation to get true anomaly `ν`
2. Calculate radius: `r = a(1-e²) / (1 + e×cos(ν))`
3. Calculate perifocal coordinates (position and velocity in orbital plane)
4. Apply rotation matrices to transform from perifocal to inertial frame:
   - Rotation by argument of periapsis (ω)
   - Rotation by inclination (i)
   - Rotation by longitude of ascending node (Ω)

**Use**: Essential for satellite propagation, mission planning, and converting between different orbital representations.

---

## Coordinate System Conventions

### Equatorial (RA/Dec)
- **Reference Frame**: Fixed relative to stars (J2000.0 epoch)
- **RA**: Right Ascension (0-24 hours, measured eastward from vernal equinox)
- **Dec**: Declination (-90° to +90°, measured from celestial equator)

### Horizontal (Alt/Az)
- **Reference Frame**: Observer-based, rotates with Earth
- **Alt**: Altitude (-90° to +90°, angle above horizon)
- **Az**: Azimuth (0-360°, measured clockwise from North)

### ECEF (Earth-Centered Earth-Fixed)
- **Reference Frame**: Rotates with Earth
- **X-axis**: Points from Earth's center through equator at prime meridian (Greenwich)
- **Y-axis**: Points from Earth's center through equator at 90°E longitude
- **Z-axis**: Points from Earth's center through North Pole
- **Units**: Meters

### ECI (Earth-Centered Inertial)
- **Reference Frame**: Fixed relative to stars (J2000.0 epoch)
- **X-axis**: Points from Earth's center toward vernal equinox
- **Y-axis**: Completes right-handed system in equatorial plane
- **Z-axis**: Points from Earth's center through North Pole
- **Units**: Meters

---

## Accuracy and Limitations

### Current Accuracy
- **Solar Position**: ~1 arcminute accuracy (suitable for most applications)
- **Lunar Position**: ~10 arcminutes accuracy (perturbation theory with limited terms)
- **Coordinate Conversions**: Sub-arcsecond accuracy (limited by floating-point precision)
- **ECEF/ECI Transformations**: ~1mm accuracy for Earth-scale coordinates

### Limitations
- **Simplified Models**: Solar and lunar positions use simplified models (not full VSOP87/ELP2000)
- **No Precession/Nutation**: ECEF/ECI transformations assume J2000.0 epoch without precession/nutation corrections
- **No Light-Time Correction**: Planetary positions (when implemented) won't include light-time correction initially
- **No Atmospheric Refraction**: Alt/Az conversions don't account for atmospheric refraction

### Future Enhancements
- Full VSOP87 implementation for planetary positions
- IAU-2000/2006 precession-nutation models for high-precision ECEF/ECI
- Atmospheric refraction corrections for Alt/Az
- Light-time and aberration corrections for planetary positions

---

## VSOP87 Planetary Position Theory

### Overview
VSOP87 (Variations Séculaires des Orbites Planétaires 1987) is a semi-analytical model developed by P. Bretagnon and G. Francou at the Bureau des Longitudes in Paris. It provides precise calculations of planetary positions in the Solar System by representing heliocentric coordinates as sums of periodic terms.

### Mathematical Foundation

#### Series Representation
Each coordinate (longitude L, latitude B, radius R) is expressed as a series of terms:

**General Form**: `A × cos(B + C × t)`

Where:
- **A**: Amplitude of the term
- **B**: Phase (constant term in radians)
- **C**: Frequency (coefficient for time variable)
- **t**: Time in Julian centuries from J2000.0

#### Series Structure
VSOP87 uses multiple series for each coordinate:

**Longitude (L)**:
```
L = (L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8
```

**Latitude (B)**:
```
B = (B0 + B1×t + B2×t² + B3×t³ + B4×t⁴) / 10^8
```

**Radius (R)**:
```
R = (R0 + R1×t + R2×t² + R3×t³ + R4×t⁴) / 10^8 (in AU)
```

Each series (L0, L1, etc.) is itself a sum of terms: `Σ(A × cos(B + C × t))`

### VSOP87 Variants

VSOP87 offers several versions for different applications:

- **VSOP87**: Heliocentric ecliptic orbital elements (equinox J2000.0)
- **VSOP87A**: Heliocentric ecliptic rectangular coordinates (equinox J2000.0)
- **VSOP87B**: Heliocentric ecliptic spherical coordinates (equinox J2000.0) - **Used in this implementation**
- **VSOP87C**: Heliocentric ecliptic rectangular coordinates (equinox of date)
- **VSOP87D**: Heliocentric ecliptic spherical coordinates (equinox of date)
- **VSOP87E**: Barycentric ecliptic rectangular coordinates (equinox J2000.0)

### Implementation Approach

**Decision**: Implement simplified VSOP87 (truncated series)

**Rationale**:
- Custom implementation demonstrates deeper understanding of the theory
- Portfolio project benefits from showing implementation skills
- Simplified version uses most significant terms (can be extended to full VSOP87)
- Full VSOP87 contains thousands of terms per planet; truncated version is more manageable

**Accuracy Expectations**:
- **Simplified VSOP87** (truncated series):
  - Inner planets (Mercury, Venus, Mars): ~1 arcminute
  - Outer planets (Jupiter, Saturn): ~1-2 arcminutes
  - Distant planets (Uranus, Neptune): ~2-5 arcminutes
- **Full VSOP87**: <1 arcsecond for all planets

### Coordinate Conversion Pipeline

1. **Heliocentric Ecliptic (L, B, R)**: Raw VSOP87 output
   - L: Ecliptic longitude (radians)
   - B: Ecliptic latitude (radians)
   - R: Radius vector (Astronomical Units)

2. **Heliocentric Equatorial**: Convert ecliptic to equatorial
   - Apply rotation using obliquity of the ecliptic
   - Time-dependent obliquity: `ε = 23.439291° - 0.0130042° × t - ...`

3. **Geocentric Equatorial**: Account for Earth's position
   - Calculate Earth's heliocentric position using VSOP87
   - Subtract Earth's position vector from planet's position
   - Optional: Apply light-time correction for high precision

4. **Final Output (RA, Dec)**: Geocentric equatorial coordinates
   - Right Ascension (hours)
   - Declination (degrees)

### Data Sources

- **Official VSOP87 Data**: Available from IMCCE (Institut de Mécanique Céleste et de Calcul des Éphémérides)
  - FTP Server: https://ftp.imcce.fr/pub/ephem/planets/vsop87/
  - Provides complete VSOP87 coefficient files for all planets
- **Truncated Coefficients**: Simplified versions for reduced computational load
  - VSOP87 Multi-Language project provides truncated versions at different precision levels
  - Truncation eliminates terms with coefficients below threshold (e.g., < 1/1000 of largest term)
- **Validation**: Compare results with JPL Horizons ephemeris service

### Coefficient Data Structure

#### Storage Strategy
**Decision**: Compile-time constants (const arrays)

**Rationale**:
- **Performance**: No runtime loading overhead, data embedded in binary
- **Simplicity**: No file I/O required, easier deployment
- **Type Safety**: Compile-time validation of data structure
- **Memory Efficiency**: Data stored as static arrays, shared across instances

**Alternative Considered**: File-based loading
- Pros: Easier to update coefficients without recompiling
- Cons: Runtime overhead, file dependency, more complex error handling
- Decision: Use compile-time constants for initial implementation; can add file loading later

#### Coefficient Format

Each VSOP87 term is stored as a `Vsop87Term` structure:

```rust
struct Vsop87Term {
    amplitude: f64,  // A: Coefficient magnitude
    phase: f64,      // B: Phase angle (radians)
    frequency: f64,  // C: Frequency term
}
```

**Example Term**:
- Amplitude: `4.40250710144` (arcseconds for L, radians for B, AU for R)
- Phase: `0.0` (radians)
- Frequency: `26087.90314157420` (for Mercury's main orbital frequency)

**Series Organization**:
- Each planet has three `Vsop87Series`: longitude (L), latitude (B), radius (R)
- Each series contains 3-6 sub-series (L0-L5, B0-B4, R0-R4)
- Each sub-series is a `Vec<Vsop87Term>` containing all terms for that series

**Data Structure Example** (Mercury - simplified):
```
PlanetVsop87Data {
    longitude: Vsop87Series {
        series_0: [term1, term2, term3, ...],  // L0: ~10-20 terms (truncated)
        series_1: [term1, term2, ...],         // L1: ~5-10 terms
        series_2: [term1, term2, ...],         // L2: ~3-5 terms
        series_3: [term1, ...],                // L3: ~1-3 terms
        series_4: Some([term1, ...]),          // L4: ~1 term
        series_5: Some([term1, ...]),          // L5: ~1 term
    },
    latitude: Vsop87Series {
        series_0: [term1, term2, ...],         // B0: ~5-10 terms
        series_1: [term1, ...],                // B1: ~1-3 terms
        series_2: [],                          // B2: minimal
        series_3: [],                          // B3: minimal
        series_4: None,                        // B4: not used
        series_5: None,                        // B5: not used
    },
    radius: Vsop87Series {
        series_0: [term1, term2, ...],         // R0: ~5-10 terms
        series_1: [term1, ...],                // R1: ~1-3 terms
        series_2: [],                          // R2: minimal
        series_3: [],                          // R3: minimal
        series_4: None,                        // R4: not used
        series_5: None,                        // R5: not used
    },
}
```

**Full vs Truncated**:
- **Full VSOP87**: Mercury has ~500-1000 terms per series (L0 alone has ~500 terms)
- **Truncated (this implementation)**: Uses ~10-20 largest terms per series
- **Accuracy Trade-off**: Truncated provides ~1 arcminute vs <1 arcsecond for full

#### Coefficient Acquisition Process

1. **Source Selection**: 
   - Primary: IMCCE official VSOP87 files
   - Alternative: VSOP87 Multi-Language project (truncated versions)

2. **Data Processing**:
   - Parse VSOP87 file format (typically space-separated: A B C)
   - Filter terms by amplitude threshold (for truncated version)
   - Organize by planet, variable (L/B/R), and series (0-5)

3. **Storage Implementation**:
   - Convert parsed data to Rust `Vsop87Term` structures
   - Organize into `Vsop87Series` and `PlanetVsop87Data`
   - Embed as const arrays in source code

4. **Validation**:
   - Verify coefficient ranges are reasonable
   - Check series structure matches expected format
   - Test with known reference values (J2000.0 epoch)

### Research Findings

**Available Rust Crates**:
- `vsop87` crate exists and provides full VSOP87 implementation
- `siderust` crate offers comprehensive astronomical computations
- `astro` crate provides planetary positioning algorithms

**Decision**: Implement custom simplified version to demonstrate:
- Understanding of VSOP87 mathematical theory
- Ability to implement complex astronomical algorithms
- Software engineering skills (data structures, algorithms, testing)

---

## References

This program implements algorithms based on:
- **Jean Meeus**: "Astronomical Algorithms" - Standard reference for astronomical calculations
- **IAU Standards**: International Astronomical Union conventions for coordinate systems
- **Kepler's Laws**: Classical orbital mechanics
- **Spherical Trigonometry**: Coordinate transformation mathematics
- **VSOP87 Theory**: P. Bretagnon and G. Francou - Planetary ephemeris theory
- **JPL Horizons**: NASA/JPL ephemeris service for validation

All calculations are verified against authoritative sources and include comprehensive test coverage.
