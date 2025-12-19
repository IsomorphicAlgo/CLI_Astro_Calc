# CLI Astro Calc - Overview

## Core Purpose

**CLI Astro Calc** is a command-line astronomy and orbital mechanics calculator built in Rust. The program provides accurate astronomical calculations for both educational and practical applications, including:

- **Time Systems**: Converting between calendar dates and astronomical time systems (Julian Date, Sidereal Time)
- **Celestial Object Positions**: Calculating the positions of the Sun, Moon, and planets in the sky
- **Rise/Set Times**: Determining when celestial objects rise and set for any location on Earth
- **Coordinate Transformations**: Converting between different astronomical coordinate systems
- **Orbital Mechanics**: Solving Kepler's equation and converting between orbital elements and state vectors
- **Earth-Centered Coordinates**: Transforming between Earth-Centered Earth-Fixed (ECEF) and Earth-Centered Inertial (ECI) coordinate systems
- **Planetary Ephemeris**: VSOP87-based planetary position calculations (heliocentric coordinates)

The program is designed for accuracy, with all calculations verified against authoritative astronomical sources and comprehensive test coverage.

---

## Top-Level Feature Summary

### 1. Time Systems (`time.rs`)
- **Julian Date Calculation**: Converts calendar dates to continuous day count since 4713 BC
- **Greenwich Mean Sidereal Time (GMST)**: Calculates Earth's rotation relative to the stars
- **Local Sidereal Time (LST)**: Converts GMST to local sidereal time for any longitude

### 2. Celestial Object Positions (`celestial.rs`, `planets.rs`)
- **Solar Position**: Calculates Sun's Right Ascension (RA) and Declination (Dec) using mean anomaly and equation of center
- **Lunar Position**: Calculates Moon's RA/Dec using perturbation theory with multiple periodic terms
- **Planetary Positions**: Calculates planet positions using VSOP87 theory (heliocentric ecliptic coordinates L, B, R)
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

The VSOP87 coordinate conversion pipeline transforms heliocentric ecliptic coordinates (L, B, R) to geocentric equatorial coordinates (RA, Dec) suitable for observations from Earth.

#### Step 1: Heliocentric Ecliptic (L, B, R) - Raw VSOP87 Output
- **L**: Ecliptic longitude (radians) - position along the ecliptic plane
- **B**: Ecliptic latitude (radians) - position above/below the ecliptic plane
- **R**: Radius vector (Astronomical Units) - distance from Sun to planet

#### Step 2: Obliquity of the Ecliptic Calculation

The obliquity of the ecliptic (ε) is the angle between Earth's equatorial plane and the ecliptic plane (Earth's orbital plane around the Sun). It changes slowly over time due to precession.

**Formula**:
```
ε = 23.4393° - 0.0000004° × d
```
Where `d` = days since J2000.0

**Higher Precision Formula** (for arcsecond accuracy):
```
ε = 23.439291° - 0.0130042° × t - 0.00000016° × t² + 0.000000503° × t³
```
Where `t` = Julian centuries from J2000.0

**Current Value**: ~23.44° (decreasing by ~0.47 arcseconds per century)

#### Step 3: Ecliptic to Equatorial Conversion

Convert from ecliptic spherical coordinates (L, B, R) to equatorial rectangular coordinates (x, y, z).

**Ecliptic Rectangular Coordinates**:
```
x_ecl = R × cos(B) × cos(L)
y_ecl = R × cos(B) × sin(L)
z_ecl = R × sin(B)
```

**Rotation to Equatorial** (rotate around X-axis by obliquity angle ε):
```
x_eq = x_ecl
y_eq = y_ecl × cos(ε) - z_ecl × sin(ε)
z_eq = y_ecl × sin(ε) + z_ecl × cos(ε)
```

This rotation accounts for the ~23.44° tilt between the ecliptic and equatorial planes.

#### Step 4: Heliocentric to Geocentric Conversion

Convert from heliocentric (Sun-centered) to geocentric (Earth-centered) coordinates by subtracting Earth's position.

**Geocentric Position Vector**:
```
[x_geo]   [x_planet]   [x_earth]
[y_geo] = [y_planet] - [y_earth]
[z_geo]   [z_planet]   [z_earth]
```

Where:
- `(x_planet, y_planet, z_planet)`: Planet's heliocentric equatorial position
- `(x_earth, y_earth, z_earth)`: Earth's heliocentric equatorial position (calculated using VSOP87)
- `(x_geo, y_geo, z_geo)`: Geocentric equatorial position

**Note**: Both planet and Earth positions must be in the same coordinate system (equatorial) before subtraction.

#### Step 5: Rectangular to RA/Dec Conversion

Convert geocentric equatorial rectangular coordinates to Right Ascension and Declination.

**Right Ascension (RA)**:
```
RA = atan2(y, x) (in radians)
RA_hours = (RA_degrees / 15.0) mod 24
```

**Declination (Dec)**:
```
r = √(x² + y² + z²)
Dec = arcsin(z / r) (in radians)
Dec_degrees = Dec_radians × (180/π)
```

Where:
- `r`: Distance from Earth to planet
- `RA`: 0-24 hours (measured eastward from vernal equinox)
- `Dec`: -90° to +90° (measured from celestial equator)

#### Complete Pipeline Summary

```
VSOP87 Output (L, B, R)
    ↓
Calculate Obliquity ε
    ↓
Ecliptic → Equatorial (rotation by ε)
    ↓
Heliocentric → Geocentric (subtract Earth position)
    ↓
Rectangular → Spherical (RA, Dec)
    ↓
Final Output (RA in hours, Dec in degrees)
```

#### Light-Time Correction (Optional)

For high-precision applications (>arcsecond accuracy), light-time correction accounts for the finite speed of light:

1. Calculate planet position at time `t`
2. Calculate distance from Earth to planet
3. Calculate light travel time: `Δt = distance / c` (where c = speed of light)
4. Recalculate planet position at time `t - Δt`
5. Use corrected position for final RA/Dec

**Impact**: Typically ~0.01 arcseconds for inner planets, negligible for outer planets.

**Current Implementation**: Light-time correction is omitted for simplicity. Can be added for high-precision applications.

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

### VSOP87 Series Evaluation

#### Time Variable Calculation

VSOP87 uses time in **Julian centuries from J2000.0** as the independent variable:

```
t = (JD - J2000.0) / 36525.0
```

Where:
- `JD`: Julian Date for the calculation
- `J2000.0`: 2451545.0 (January 1, 2000, 12:00:00 TT)
- `36525.0`: Number of days in a Julian century

**Example**: For January 1, 2024 (JD ≈ 2460311.0):
```
t = (2460311.0 - 2451545.0) / 36525.0 ≈ 0.2400 centuries
```

#### Series Evaluation Algorithm

Each VSOP87 coordinate (L, B, or R) is calculated by evaluating multiple series and combining them with time powers:

**Step 1: Evaluate Individual Terms**

For each term in a series, calculate:
```
term_value = A × cos(B + C × t)
```

Where:
- `A`: Amplitude (coefficient)
- `B`: Phase (constant term in radians)
- `C`: Frequency (coefficient for time variable)
- `t`: Time in Julian centuries from J2000.0

**Step 2: Sum Terms in Each Series**

For each sub-series (L0, L1, L2, etc.), sum all terms:
```
L0 = Σ(A_i × cos(B_i + C_i × t))  for all terms i in series_0
L1 = Σ(A_i × cos(B_i + C_i × t))  for all terms i in series_1
...
```

**Step 3: Combine Series with Time Powers**

Combine all sub-series using polynomial in time:
```
L = (L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8
B = (B0 + B1×t + B2×t² + B3×t³ + B4×t⁴) / 10^8
R = (R0 + R1×t + R2×t² + R3×t³ + R4×t⁴) / 10^8
```

**Note**: The division by 10^8 is because VSOP87 coefficients are stored in units where the final result needs scaling.

#### Implementation Details

**Longitude (L) Calculation**:
- Uses 6 sub-series: L0, L1, L2, L3, L4, L5
- Formula: `L = (L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8`
- Result in radians (converted from arcseconds in original VSOP87)
- Normalized to [0, 2π) range

**Latitude (B) Calculation**:
- Uses 5 sub-series: B0, B1, B2, B3, B4
- Formula: `B = (B0 + B1×t + B2×t² + B3×t³ + B4×t⁴) / 10^8`
- Result in radians
- Typically small values (planets stay near ecliptic plane)

**Radius (R) Calculation**:
- Uses 5 sub-series: R0, R1, R2, R3, R4
- Formula: `R = (R0 + R1×t + R2×t² + R3×t³ + R4×t⁴) / 10^8`
- Result in Astronomical Units (AU)
- Represents distance from Sun to planet

#### Performance Considerations

**Optimization Strategies**:
1. **Pre-compute time powers**: Calculate t², t³, t⁴, t⁵ once and reuse
2. **Vectorized operations**: Evaluate multiple terms in parallel (SIMD if available)
3. **Term filtering**: For truncated versions, only evaluate significant terms
4. **Caching**: Cache intermediate results if calculating multiple planets at same epoch

**Computational Complexity**:
- For truncated VSOP87 (N terms per series):
  - Time per series: O(N) for term evaluation
  - Total time: O(N × M) where M = number of sub-series (typically 3-6)
  - For full VSOP87: N can be 500-1000 terms per series

**Typical Performance**:
- Truncated VSOP87: < 1ms per planet per epoch
- Full VSOP87: ~5-10ms per planet per epoch
- Optimized implementations can achieve < 0.1ms for truncated versions

#### Validation Methodology

**Test Strategy**:

1. **Unit Tests for VSOP87 Series Evaluation**:
   - Test with known coefficient sets and verify trigonometric calculations
   - Test edge cases: t = 0 (J2000.0 epoch), large t values (±20 centuries), negative t values
   - Verify series evaluation is deterministic (same input produces same output)
   - Test empty series and series with single terms

2. **Unit Tests for Planet Coordinates**:
   - Test heliocentric ecliptic coordinates (L, B, R) at J2000.0 epoch
   - Verify coordinate ranges:
     - Longitude L: [0, 2π) radians
     - Latitude B: [-π/2, π/2] radians (planets stay near ecliptic)
     - Radius R: Positive values in Astronomical Units
   - Test multiple epochs (past and future) to verify numerical stability

3. **Integration Tests**:
   - Test multiple planets at the same epoch (verifies consistency)
   - Test planets at different epochs (past and future dates)
   - Verify that positions change over time (planets move in their orbits)
   - Test coordinate normalization (longitude wraps correctly)

4. **Reference Sources for Validation**:
   - **JPL Horizons**: NASA/JPL ephemeris service provides authoritative planetary positions
     - URL: https://ssd.jpl.nasa.gov/horizons/
     - Provides high-precision positions for validation
   - **Jean Meeus "Astronomical Algorithms"**: Standard reference for astronomical calculations
     - Contains VSOP87 examples and validation data
   - **IMCCE VSOP87**: Official VSOP87 data source
     - Provides reference implementations and test cases

5. **Accuracy Requirements**:
   - **Simplified VSOP87 (truncated series)**:
     - Inner planets (Mercury, Venus, Mars): ~1 arcminute accuracy
     - Outer planets (Jupiter, Saturn): ~1-2 arcminutes accuracy
     - Distant planets (Uranus, Neptune): ~2-5 arcminutes accuracy
   - **Full VSOP87**:
     - All planets: < 1 arcsecond accuracy
   - **Validation Approach**:
     - Compare calculated positions with JPL Horizons at J2000.0 epoch
     - Verify differences are within expected accuracy ranges
     - Test at multiple epochs to ensure consistency

6. **Performance Benchmarks**:
   - Target: < 10ms per planet calculation (for truncated VSOP87)
   - Actual: < 1ms per planet for truncated VSOP87 (measured)
   - Series evaluation: < 100 microseconds per series evaluation
   - Performance tests run 1000+ iterations to get reliable averages

7. **Error Handling and Edge Cases**:
   - Invalid Julian Date (NaN, infinity) → Error with descriptive message
   - Extreme dates (> ±20 centuries from J2000.0) → Warning logged, calculation proceeds
   - Missing planet data → Error with planet name
   - Coordinate validation → Ensures all coordinates are in valid ranges

8. **Test Coverage**:
   - Unit tests: VSOP87 term evaluation, series evaluation, coordinate calculations
   - Integration tests: Multiple planets, multiple epochs, coordinate ranges
   - Performance tests: Calculation speed, series evaluation speed
   - Validation tests: Input validation, error handling, edge cases
   - Current test count: 15+ comprehensive tests covering all major functionality

**Validation Process**:
1. Implement VSOP87 evaluation functions
2. Write unit tests for individual components
3. Write integration tests for complete workflows
4. Compare results with reference sources (JPL Horizons)
5. Verify accuracy within expected ranges
6. Performance testing to ensure acceptable speed
7. Document validation results and any limitations

**Future Enhancements**:
- Automated comparison with JPL Horizons API
- Regression tests with known-good reference values
- Continuous integration with validation checks
- Extended accuracy testing for full VSOP87 implementation

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

## CLI Usage for Planetary Positions

### Position Command

The `position` command calculates the Right Ascension (RA) and Declination (Dec) of celestial objects, including planets.

**Basic Syntax**:
```bash
cargo run -- position --object <OBJECT> --date <YYYY-MM-DD>
```

**Supported Objects**:
- **Sun**: `--object "sun"`
- **Moon**: `--object "moon"`
- **Planets**: `--object "mercury"`, `"venus"`, `"mars"`, `"jupiter"`, `"saturn"`, `"uranus"`, `"neptune"`

**Example Usage**:
```bash
# Calculate Jupiter's position on January 1, 2000
cargo run -- position --object "jupiter" --date "2000-01-01"

# Calculate Mars position on a specific date
cargo run -- position --object "mars" --date "2024-12-25"
```

**Output Format**:
```
RA:  HH:MM:SS
Dec: ±DD°MM'SS"
```

Where:
- **RA**: Right Ascension in hours, minutes, seconds (0-24 hours)
- **Dec**: Declination in degrees, arcminutes, arcseconds (-90° to +90°)

### Calculation Pipeline

When calculating a planet's position, the CLI follows this pipeline:

1. **Parse Input**: Date string → Julian Date
2. **VSOP87 Evaluation**: Calculate heliocentric ecliptic coordinates (L, B, R)
3. **Coordinate Conversion**: 
   - Ecliptic → Equatorial (rotation by obliquity)
   - Heliocentric → Geocentric (subtract Earth's position)
   - Rectangular → Spherical (RA, Dec)
4. **Format Output**: Convert to hours/degrees format

### Error Handling

**Common Errors**:
- **"Unknown object"**: Planet name not recognized (check spelling, case-insensitive)
- **"Earth VSOP87 data not available"**: Earth's position data needed for geocentric conversion (may occur if Earth data is placeholder)
- **"Invalid date"**: Date format must be YYYY-MM-DD

**Verbose Logging**:
Use `--verbose` flag to see detailed calculation steps:
```bash
cargo run -- --verbose position --object "jupiter" --date "2000-01-01"
```

This shows:
- VSOP87 series evaluation values
- Heliocentric ecliptic coordinates
- Obliquity calculation
- Coordinate conversion steps
- Final RA/Dec values

### Accuracy Considerations

**For Planets**:
- **Truncated VSOP87**: ~1-5 arcminutes accuracy (current implementation)
- **Full VSOP87**: <1 arcsecond accuracy (if full coefficients are used)
- **Light-time correction**: Not applied (adds ~0.01 arcseconds for inner planets)

**For Sun/Moon**:
- **Sun**: ~1 arcminute accuracy (simplified model)
- **Moon**: ~10 arcminutes accuracy (perturbation theory with limited terms)

### Rise/Set Times

The `rise-set` command supports planets, but planet rise/set calculations are not yet fully implemented. The command will return an error indicating that planet rise/set times are pending.

**Future Enhancement**: Planet rise/set times will use the same algorithm as Sun/Moon, but with planet positions calculated via VSOP87.

---

## Error Handling and Logging Patterns

### Logging Strategy

The VSOP87 implementation uses a multi-level logging system to provide visibility into calculations while maintaining performance:

**Debug Level** (`--verbose` flag):
- VSOP87 series evaluation details (term counts, intermediate values)
- Obliquity calculations
- Coordinate conversion steps (ecliptic → equatorial, heliocentric → geocentric)
- Rectangular coordinate values at each step

**Info Level** (default):
- Planet calculation start/completion
- Time in Julian centuries from J2000.0
- Heliocentric ecliptic coordinates (L, B, R)
- Earth's heliocentric position
- Final geocentric RA/Dec coordinates

**Warn Level**:
- Extreme dates (>±20 centuries from J2000.0)
- Out-of-range Julian Dates
- Potential accuracy issues
- Placeholder data detection
- Unusual coordinate values

**Error Level**:
- Invalid input validation failures
- Calculation errors
- Missing data errors

### Error Handling Patterns

#### Input Validation

**Julian Date Validation**:
```rust
// Check for NaN and infinity
if julian_date.is_nan() {
    return Err("Julian Date cannot be NaN");
}
if julian_date.is_infinite() {
    return Err("Julian Date cannot be infinite");
}

// Check reasonable range
if julian_date < MIN_JD || julian_date > MAX_JD {
    warn!("Julian Date outside recommended range");
}
```

**Time Argument Validation**:
- Valid range: ~2000 BC to 3000 AD (JD ~1,000,000 to ~3,000,000)
- VSOP87 optimal range: ±20 centuries from J2000.0
- Warnings issued for extreme dates, but calculation proceeds

#### Data Availability Checks

**VSOP87 Coefficient Data**:
```rust
let vsop87_data = get_planet_vsop87_data(planet)
    .ok_or_else(|| {
        error!("VSOP87 data not available for {}", planet.name());
        AstroError::InvalidCoordinate(
            format!("VSOP87 data not available for {}. This planet may not be fully implemented yet.", planet.name())
        )
    })?;
```

**Placeholder Data Detection**:
- Checks if Earth's VSOP87 data is placeholder (empty series)
- Warns user but allows calculation to proceed
- Helps identify when full VSOP87 data needs to be added

#### Calculation Result Validation

**Coordinate Range Validation**:
- Longitude: Normalized to [0, 2π) radians
- Latitude: Validated to [-π/2, π/2] radians (warns if outside)
- Radius: Must be positive (error if ≤ 0)
- RA: Validated to [0, 24) hours
- Dec: Validated to [-90°, +90°] degrees

**NaN Detection**:
- Checks for NaN values in all calculated coordinates
- Returns error with descriptive message
- Helps identify calculation issues early

#### Error Recovery Approaches

**Graceful Degradation**:
- Extreme dates: Warn but continue calculation
- Placeholder data: Warn but attempt calculation
- Out-of-range coordinates: Warn but return result

**Fail-Fast for Critical Errors**:
- Invalid Julian Date (NaN/infinity): Return error immediately
- Missing required data (Earth position): Return error immediately
- Negative radius: Return error immediately

**Error Messages**:
- Descriptive: Explain what went wrong
- Actionable: Suggest how to fix the issue
- Contextual: Include relevant values (planet name, Julian Date, etc.)

### Example Error Scenarios

**Scenario 1: Missing Planet Data**
```
Error: VSOP87 data not available for Jupiter. This planet may not be fully implemented yet.
```
**Recovery**: Ensure planet's VSOP87 coefficients are added to `get_jupiter_vsop87_data()`

**Scenario 2: Missing Earth Data**
```
Error: Earth VSOP87 data not available for geocentric conversion. 
Earth's position is required to convert from heliocentric to geocentric coordinates.
```
**Recovery**: Implement Earth's VSOP87 coefficients

**Scenario 3: Extreme Date**
```
Warning: Julian Date 3000000.0 is 15.04 centuries from J2000.0. 
VSOP87 accuracy may degrade for extreme dates (>±20 centuries).
```
**Recovery**: Calculation proceeds, but user is warned about potential accuracy issues

**Scenario 4: Invalid Input**
```
Error: Julian Date cannot be NaN (Not a Number). Please provide a valid date.
```
**Recovery**: User must provide a valid Julian Date

### Best Practices

1. **Validate Early**: Check inputs at function entry
2. **Log Context**: Include relevant values in log messages
3. **Warn for Degradation**: Alert users to potential accuracy issues
4. **Fail Fast for Critical Errors**: Don't proceed with invalid data
5. **Provide Actionable Errors**: Tell users how to fix the problem
6. **Use Appropriate Log Levels**: Debug for details, Info for operations, Warn for issues, Error for failures

---

## Testing and Validation Methodology

### Test Coverage Strategy

The CLI Astro Calc project uses a comprehensive multi-level testing approach to ensure accuracy and reliability:

**Unit Tests** (71 tests):
- Individual function testing with known reference values
- Edge case validation (poles, equator, origin, large coordinates)
- Input validation (NaN, infinity, range checks)
- Round-trip accuracy tests (verify transformations are reversible)
- Performance benchmarks (ensure calculations meet speed targets)

**Integration Tests**:
- End-to-end workflow testing (CLI commands with various inputs)
- Coordinate system conversion pipelines
- Multiple planet calculations at same/different epochs
- Cross-module functionality (time → coordinates → celestial objects)

**Validation Tests**:
- Comparison with authoritative sources (JPL Horizons, Meeus algorithms)
- Known coordinate pair verification (Greenwich meridian, poles, equator)
- Numerical stability testing (large coordinates, extreme dates)
- Accuracy verification (within acceptable tolerances)

**Doctests** (5 tests):
- Code examples in documentation are verified to compile and run
- Ensures documentation stays synchronized with implementation

### Accuracy Validation

**ECEF/ECI Transformations**:
- Round-trip accuracy: < 1mm for Earth-scale coordinates (~6,378 km)
- Verified at multiple GMST values (0, 6, 12, 18, 23.5 hours)
- Tested with known coordinate pairs (Greenwich meridian, poles, equator)
- Numerical stability verified for geosynchronous orbit distances (~42,164 km)

**VSOP87 Planetary Positions**:
- Mercury position verified at J2000.0 epoch (reference: JPL Horizons)
- Coordinate range validation (longitude 0-2π, latitude ±π/2, radius > 0)
- Time calculation accuracy (Julian centuries from J2000.0)
- Performance targets: < 1ms per planet calculation (actual: < 1ms achieved)

**Coordinate Conversions**:
- RA/Dec ↔ Alt/Az round-trip accuracy verified
- Edge cases tested (zenith, horizon, poles, equator)
- Coordinate normalization validated

### Performance Benchmarks

**Target vs Actual Performance**:
- Planet calculations: Target < 10ms, Actual < 1ms ✅
- VSOP87 series evaluation: Target < 100μs, Actual < 100μs ✅
- ECEF/ECI transformations: < 1μs per transformation ✅
- Coordinate conversions: < 1μs per conversion ✅

**Optimization Techniques**:
- Pre-computed time powers for VSOP87 (avoid repeated calculations)
- Efficient rotation matrix operations (direct matrix multiplication)
- Compile-time constant storage for VSOP87 coefficients (no runtime loading)

### Test Execution

**Running Tests**:
```bash
# Run all tests
cargo test

# Run only unit tests
cargo test --lib

# Run only doctests
cargo test --doc

# Run with verbose output
cargo test -- --nocapture

# Run performance tests
cargo test --release test_performance
```

**Test Organization**:
- Tests are co-located with source code (in `#[cfg(test)]` modules)
- Each module has its own test module (`mod tests { ... }`)
- Integration tests verify complete workflows
- Performance tests ensure no regressions

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
