# CLI Astro Calc

A command-line astronomy and orbital mechanics calculator built in Rust.

## Version 0.1

**Features:**
- Julian Date and Sidereal Time calculations
- Solar and lunar position calculations (RA/Dec)
- Rise/set times for Sun and Moon
- Coordinate conversions (RA/Dec ↔ Alt/Az)
- Kepler's equation solver
- Orbital period calculations
- Orbital elements to state vectors

### Phase 1 - Basic Structure ✅ (Completed)
- ✅ Project structure with proper Rust organization
- ✅ CLI interface with subcommands
- ✅ Logging and error handling systems
- ✅ Module structure and documentation

### Phase 2 - Core Astronomy Calculations ✅ (Completed)
- ✅ **Julian Date calculation** - Foundation for all astronomical calculations
- ✅ **Greenwich Mean Sidereal Time (GMST)** - Earth's rotation relative to stars
- ✅ **Local Sidereal Time (LST)** - Sidereal time at any longitude
- ✅ **Solar Position** - Calculate Sun's RA/Dec coordinates for any date/time
- ✅ **Rise/Set Times** - Calculate sunrise/sunset times for any location
- ✅ **Comprehensive Testing** - All calculations verified with astronomical accuracy

### Phase 3 - Coordinate Conversions ✅ (Completed)
- ✅ **RA/Dec → Alt/Az** - Convert equatorial to horizontal coordinates
- ✅ **Alt/Az → RA/Dec** - Convert horizontal to equatorial coordinates
- ✅ **CLI Integration** - Full command-line interface for conversions
- ✅ **Real-time LST** - Uses current time for accurate conversions
- ✅ **Comprehensive Testing** - Edge cases and round-trip accuracy verified

### Phase 4 - Lunar & Orbital Mechanics ✅ (Completed)
- ✅ **Lunar Position** - Calculate Moon's RA/Dec using perturbation theory
- ✅ **Lunar Rise/Set** - Calculate moonrise and moonset times
- ✅ **Kepler's Equation** - Solve for true anomaly from mean anomaly
- ✅ **Orbital Period** - Calculate period using Kepler's Third Law
- ✅ **Elements → State Vectors** - Convert orbital elements to position/velocity
- ✅ **Comprehensive Testing** - 29 tests covering all functions

### Still To Implement
- ECEF ↔ ECI transformations
- Planet positions (using VSOP87 or similar)
- More advanced orbital propagation

## Project Structure

```
CLI_Astro_Calc/
├── Cargo.toml
├── src/
│   ├── main.rs         # CLI entry point
│   ├── lib.rs          # Library and error types
│   ├── coordinates.rs  # Coordinate conversions
│   ├── celestial.rs    # Sun/Moon calculations
│   ├── time.rs         # Time systems
│   └── orbital.rs      # Orbital mechanics
└── readme.md
```

## Dependencies

- clap: CLI parsing
- chrono: Date/time handling
- thiserror/anyhow: Error handling
- log/env_logger: Logging

## Usage

```bash
# Build
cargo build

### Run
```bash
# Show help
cargo run -- --help

# Working examples with real calculations:

# Calculate sunrise/sunset times for New York
cargo run -- rise-set --object "sun" --latitude 40.7128 --longitude=-74.0060
# Output: Rise time: 2025-10-26 11:24:17 UTC, Set time: 2025-10-26 22:04:26 UTC

# Calculate solar position for a specific date
cargo run -- position --object "sun" --date "2024-06-21"
# Output: RA: 06:02:52, Dec: +23:26:02 (summer solstice)

# Calculate Julian Date and sidereal time
cargo run -- time --date "2024-01-01" --time "12:00:00"
cargo run -- convert --from ra-dec --to alt-az --coords "12.5,45.0"

# Test
cargo test
```

## Formula Reference

### Time Systems
- [Julian Date](https://en.wikipedia.org/wiki/Julian_day) - Continuous day count since 4713 BC
- [Sidereal Time](https://en.wikipedia.org/wiki/Sidereal_time) - Earth's rotation relative to stars

### Coordinate Conversions
- **RA/Dec to Alt/Az**: Convert from equatorial coordinates (Right Ascension/Declination) to horizontal coordinates (Altitude/Azimuth) for a given observer location and time
- **Alt/Az to RA/Dec**: Convert from horizontal coordinates to equatorial coordinates
- **ECEF/ECI**: Convert between Earth-Centered Earth-Fixed and Earth-Centered Inertial coordinate systems
### Celestial Coordinates
- [Equatorial Coordinate System](https://en.wikipedia.org/wiki/Equatorial_coordinate_system) - Right Ascension / Declination
- [Horizontal Coordinate System](https://en.wikipedia.org/wiki/Horizontal_coordinate_system) - Altitude / Azimuth
- [Spherical Trigonometry](https://en.wikipedia.org/wiki/Spherical_trigonometry) - Coordinate transformations

### Solar Calculations
- [Position of the Sun](https://en.wikipedia.org/wiki/Position_of_the_Sun) - Solar coordinates
- [Sunrise Equation](https://en.wikipedia.org/wiki/Sunrise_equation) - Rise/set times
- [Equation of Time](https://en.wikipedia.org/wiki/Equation_of_time) - Solar corrections

### Lunar Calculations
- [Lunar Theory](https://en.wikipedia.org/wiki/Lunar_theory) - Moon's complex orbit
- [Jean Meeus Algorithms](https://en.wikipedia.org/wiki/Jean_Meeus) - Astronomical calculations

### Celestial Calculations
- **Rise/Set Times**: Calculate when celestial objects rise and set for a given observer location
- **Solar/Lunar Positions**: Compute the current position of the Sun or Moon in the sky

### Time Calculations
- **Julian Dates**: Convert calendar dates to Julian Date system used in astronomy
- **Sidereal Time**: Calculate sidereal time (star time) for astronomical observations

### Orbital Mechanics (References)
- [Kepler's Laws](https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion) - Orbital period and motion
- [Kepler's Equation](https://en.wikipedia.org/wiki/Kepler%27s_equation) - Mean to true anomaly
- [Orbital Elements](https://en.wikipedia.org/wiki/Orbital_elements) - Classical six elements
- [Orbital State Vectors](https://en.wikipedia.org/wiki/Orbital_state_vectors) - Position and velocity

## Development Notes

- All calculations are currently placeholder implementations
- Error handling is set up with custom error types
- Logging is configured for debugging and monitoring
- The project follows Rust best practices with proper module organization
- Each calculation module is documented with function descriptions

## Phase 2 Achievements

Phase 2 has been successfully completed with the following implementations:

### Implemented Calculations
1. **Julian Date** - Converts calendar dates to continuous day count (astronomical standard)
2. **Sidereal Time** - GMST and LST calculations for stellar observations
3. **Solar Position** - Accurate RA/Dec coordinates using astronomical algorithms
4. **Rise/Set Times** - Calculates sunrise/sunset for any location on Earth

### Astronomical Accuracy
- **Seasonal Variation**: Correctly models Sun's annual motion
- **Solstice Values**: Winter (-23.4°) and Summer (+23.4°) declinations verified
- **Polar Regions**: Handles midnight sun and polar night conditions
- **Refraction**: Accounts for atmospheric light bending (-0.833°)

### Testing
All calculations have comprehensive unit tests verifying astronomical accuracy:
- Julian Date: J2000.0 epoch verification
- Sidereal Time: GMST/LST calculations  
- Solar Position: Solstice and equinox validation
- Rise/Set Times: Range and polar condition checks

## Phase 3 Achievements

Phase 3 successfully implements coordinate system transformations:

### Coordinate Conversions
1. **RA/Dec → Alt/Az** - Transforms equatorial coordinates (celestial sphere) to horizontal coordinates (observer's local sky)
   - Uses spherical trigonometry for accurate transformations
   - Accounts for observer latitude and Local Sidereal Time
   - Essential for telescope pointing and sky observation

2. **Alt/Az → RA/Dec** - Inverse transformation from local to celestial coordinates
   - Enables tracking celestial objects from local observations
   - Critical for converting telescope positions to star catalog coordinates

### Key Features
- **Real-time Calculations**: Uses current time and LST for accurate conversions
- **Observer Location**: Defaults to New York (40.7128°N, 74.0060°W) for demonstrations
- **Formatted Output**: Displays coordinates in both decimal degrees/hours and HMS/DMS format
- **Edge Case Handling**: Properly handles horizon, zenith, and pole singularities

### Testing & Validation
- Comprehensive unit tests for both conversion directions
- Zenith and horizon position tests
- Round-trip accuracy verification
- Coordinate range validation

## Phase 4 Achievements

Phase 4 successfully implements lunar calculations and orbital mechanics:

### Lunar Calculations
1. **Lunar Position** - Accurate Moon position using Jean Meeus's algorithms
   - Accounts for major perturbations from Sun and Earth
   - Includes 10+ perturbation terms for longitude and latitude
   - Achieves ~0.1° accuracy for most purposes
   
2. **Lunar Rise/Set Times** - Calculate moonrise and moonset
   - Similar algorithm to solar rise/set but accounts for Moon's faster motion (~13°/day)
   - Handles different lunar angular diameter (0.25° vs Sun's 0.267°)
   - Works at all latitudes (polar regions handled correctly)

### Orbital Mechanics
1. **Kepler's Equation Solver** - Converts mean anomaly to true anomaly
   - Uses Newton-Raphson iteration for accurate convergence
   - Handles circular to highly eccentric orbits (e < 1)
   - Converges to 10⁻¹⁰ precision in typically < 5 iterations

2. **Orbital Period** - Kepler's Third Law implementation
   - T = 2π√(a³/μ) for any two-body system
   - Verified with Earth's orbit (~1 year) and ISS orbit (~90 minutes)
   
3. **Elements to State Vectors** - Full orbital element conversion
   - Converts 6 classical orbital elements to position/velocity vectors
   - Implements proper rotation matrices (perifocal to inertial frame)
   - Conserves orbital energy (verified in tests)
   - Handles all orbital geometries (inclined, eccentric, etc.)

### Key Features
- **Perturbation Theory**: Lunar calculations include solar perturbations
- **Energy Conservation**: Orbital mechanics preserve physical laws
- **Robust Numerics**: Newton-Raphson with proper convergence criteria
- **Comprehensive Testing**: 29 unit tests covering all edge cases

### Physical Constants
- Earth's GM: 398,600 km³/s²
- Sun's GM: 1.32712440018×10¹¹ km³/s²
- J2000.0 epoch: JD 2451545.0 (Jan 1, 2000 12:00 UTC)

## Next Phase
Future phases will implement:
- ECEF ↔ ECI transformations (Earth-fixed to inertial frames)
- Planet positions (using VSOP87 or simplified algorithms)
- Orbital propagation (SGP4 for satellites)
- More coordinate systems (Galactic, Ecliptic)

---

**Important**: Remember to commit and push to GitHub after each phase completion!
