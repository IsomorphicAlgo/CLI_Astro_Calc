# CLI Astro Calc

A command-line astronomy and orbital mechanics calculator built in Rust.

## Project Overview

This project is part of a two-component portfolio system designed to showcase skills for the **space industry**. The goal is to demonstrate proficiency in astronomical calculations, orbital mechanics, and space systems engineering through practical, production-ready tools.

### Project Structure

**Part 1: CLI Tool for Astronomical Calculations** ✅ (Nearly Complete)
- A comprehensive command-line Swiss Army knife for space calculations
- Features coordinate system conversions, celestial object positions, rise/set times, and orbital mechanics
- Built with Rust for performance and reliability
- Fully tested with comprehensive test coverage

**Part 2: Space Weather Web Service** ⏳ (Coming Next)
- REST API service for fetching, caching, and serving space weather data
- Data sources: NOAA Space Weather API and similar sources
- Features: local caching, historical data storage, rate limiting, authentication
- Deployment: Self-hosted on personal server rack

### Server Infrastructure

**Hardware Specifications:**
- **CPU**: 2x 8-core 8-thread Xeon processors
- **Memory**: 32GB DDR4 ECC (Error-Correcting Code)
- **Storage**: SAS3 12-drive backplane
- **Network**: 4x 10G RJ45 ports
- **Management**: IPMI (Intelligent Platform Management Interface)
- **Power**: Redundant 800W PSUs
- **Expansion**: Room for multiple GPUs if needed

**Operating System Requirements:**
- Need to select or create an operating system suitable for hosting:
  - Rust CLI application (Part 1)
  - REST API service (Part 2 - Space Weather)
  - Database (SQLite or PostgreSQL for historical data)
  - Web server capabilities
  - Container support (optional, for deployment flexibility)
- Considerations: Linux distribution (Ubuntu Server, Debian, or custom), containerization (Docker/Kubernetes), or bare-metal deployment

### Purpose

This project serves as a portfolio demonstration of:
- **Astronomical Calculations**: Accurate implementation of standard algorithms (Julian Date, sidereal time, coordinate transformations)
- **Orbital Mechanics**: Kepler's equation, orbital elements, state vector conversions
- **Software Engineering**: Clean architecture, comprehensive testing, error handling, logging
- **Space Systems Knowledge**: Understanding of coordinate systems (ECEF/ECI, RA/Dec, Alt/Az), time systems, and celestial mechanics

All calculations are verified against authoritative sources and include comprehensive documentation.

**📖 For detailed mathematical equations and coordinate system explanations, see [OVERVIEW.md](OVERVIEW.md)**

---

## Version 0.1 - CLI Tool Status

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

## Session 2 - ECEF/ECI Transformations & Planet Positions

### Overview
Session 2 focuses on implementing two major features:
1. **ECEF ↔ ECI Transformations** - Earth-Centered Earth-Fixed to Earth-Centered Inertial coordinate system conversions
2. **Planet Positions** - Calculate positions of major planets using VSOP87 theory

### Iterative Implementation Plan

#### Part A: ECEF ↔ ECI Transformations

**Step A1: Research & Design** ✅ (Completed)
- ✅ Research ECEF/ECI transformation mathematics
  - ✅ Understand rotation matrix based on GMST
  - ✅ Review IAU-76/FK5 precession-nutation models (documented for future enhancement)
  - ✅ Document coordinate system conventions (J2000.0 epoch)
- ✅ Design function signatures in `coordinates.rs`
  - ✅ `ecef_to_eci(ecef: Ecef, gmst: f64) -> Result<Eci>` - Fully implemented with validation
  - ✅ `eci_to_ecef(eci: Eci, gmst: f64) -> Result<Ecef>` - Fully implemented with validation
  - ✅ Helper functions: `rotation_matrix_z()` and `apply_rotation_matrix()`
- ✅ Define test cases and expected results
  - ✅ Test structure created for known ECEF/ECI pairs validation
  - ✅ Round-trip accuracy tests defined
  - ✅ Edge cases defined (poles, equator, prime meridian, origin, large coordinates)
  
**Implementation Details:**
- Comprehensive documentation added to all functions
- Input validation (NaN, infinity checks)
- GMST normalization (0-24 hour range)
- Coordinate system conventions documented (J2000.0 epoch, units in meters)
- 15 test cases structured (to be implemented in Step A3)

**Step A2: Core Implementation** ✅ (Completed)
- ✅ Implement rotation matrix calculation
  - ✅ Created `rotation_matrix_z(angle_rad: f64) -> [[f64; 3]; 3]` helper
  - ✅ Uses GMST in radians for rotation angle
  - ✅ Logs rotation matrix construction (debug level)
- ✅ Implement `ecef_to_eci` function
  - ✅ Applies rotation matrix to ECEF coordinates
  - ✅ Input validation (NaN, infinity checks)
  - ✅ Logging at info level for transformations (input/output coordinates, GMST)
- ✅ Implement `eci_to_ecef` function
  - ✅ Applies rotation matrix to ECI coordinates (inverse transformation)
  - ✅ Numerical stability ensured (normalized GMST, proper matrix operations)
  - ✅ Logging at info level for transformations (input/output coordinates, GMST)
  
**Implementation Details:**
- Rotation matrix helper function with debug logging (angle, matrix values)
- Both transformation functions log at info level:
  - Input coordinates and GMST (with normalization)
  - Output coordinates
- All logging uses appropriate precision for coordinate values

**Step A3: Testing & Validation** ✅ (Completed)
- ✅ Write unit tests for rotation matrix
  - ✅ Test 0°, 90°, 180°, 270° rotations (all pass)
  - ✅ Verify matrix properties (orthogonal, determinant = 1) for multiple angles
- ✅ Write unit tests for `ecef_to_eci`
  - ✅ Test known coordinate pairs (Greenwich meridian, poles, equator)
  - ✅ Test at Greenwich meridian (x-axis in ECEF) at GMST=0 and GMST=6
  - ✅ Test at different GMST values (0, 6, 12, 18, 23.5 hours)
  - ✅ Test GMST normalization (25→1, -1→23)
  - ✅ Test edge cases (origin, large coordinates, invalid inputs)
- ✅ Write unit tests for `eci_to_ecef`
  - ✅ Test round-trip accuracy (ECEF → ECI → ECEF)
  - ✅ Verify accuracy within acceptable tolerance (< 1mm for Earth radius scale)
  - ✅ Test round-trip at multiple GMST values
- ✅ Write integration tests
  - ✅ Test with various coordinate configurations
  - ✅ Test at multiple time epochs (different GMST values)
  - ✅ Verify numerical stability for large coordinates (geosynchronous orbit)
  
**Test Results:**
- 15 comprehensive test cases implemented
- All tests verify accuracy within 1mm tolerance for Earth-scale coordinates
- Error handling tests for NaN, infinity, and invalid inputs
- Round-trip tests confirm transformation accuracy

**Step A4: CLI Integration** ✅ (Completed)
- ✅ Add `ecef-eci` conversion option to `convert` command
  - ✅ Extended `Commands::Convert` enum with `--gmst` option
  - ✅ Added coordinate parsing for "x,y,z" format
  - ✅ Added `--gmst` option (optional, auto-calculates from current time if not provided)
- ✅ Update help text and documentation
- ✅ Added example usage to readme
- ⏳ Test CLI with various inputs (ready for testing)

**Step A5: Logging & Error Handling** ✅ (Completed)
- ✅ Add comprehensive logging
  - ✅ Debug: Rotation matrix values, intermediate calculations (rotation angle, matrix application)
  - ✅ Info: Transformation operations, coordinate systems, distance calculations
  - ✅ Warn: Potential numerical issues, large coordinate values, distance changes
- ✅ Enhance error handling
  - ✅ Validate input coordinates (reasonable Earth-centered ranges, warn for > 1,000,000 km)
  - ✅ Handle edge cases (coordinates at origin, very large values with warnings)
  - ✅ Provide meaningful error messages (descriptive messages for NaN, infinity, and invalid inputs)
  
**Implementation Details:**
- Debug logging added for rotation angle calculations and matrix application steps
- Warning logs for coordinates exceeding 1,000,000 km (may indicate input errors)
- Warning logs for coordinates very close to origin (< 1 km, may be intentional for testing)
- Warning logs for distance changes during transformation (potential numerical precision issues)
- Enhanced error messages with context and suggestions
- Distance calculations logged at info level for both input and output coordinates

---

---

#### Part B: Planet Positions (VSOP87)

**Step B1: Research & Design** ✅ (Completed)
- ✅ Research VSOP87 theory
  - ✅ Understand VSOP87 series representation: A × cos(B + C × t) where t is Julian centuries from J2000.0
  - ✅ Review available VSOP87 data sources: Official VSOP87 data files available from IMCCE
  - ✅ Decide on implementation approach:
    - **Decision: Implement simplified VSOP87 (truncated series)**
    - Rationale: For portfolio demonstration, implementing our own (even simplified) shows deeper understanding
    - Alternative `vsop87` crate exists but custom implementation better demonstrates skills
    - Will use truncated series (fewer terms) for initial implementation, can be extended later
    - Full VSOP87 contains thousands of terms per planet; simplified version uses most significant terms
- ✅ Design data structures
  - ✅ `Planet` enum (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune)
  - ✅ `Vsop87Term` struct for individual series terms (amplitude, phase, frequency)
  - ✅ `Vsop87Series` struct for series L0-L5, B0-B4, R0-R4
  - ✅ `PlanetVsop87Data` struct containing longitude, latitude, and radius series
  - ✅ `HeliocentricEcliptic` struct for intermediate coordinates (L, B, R)
  - ✅ Coefficient storage: Will use const arrays for compile-time data (efficient, no runtime loading)
- ✅ Design function signatures in `planets.rs`
  - ✅ `calculate_planet_position(planet: Planet, jd: f64) -> Result<RaDec>` - Main public API
  - ✅ `evaluate_vsop87_series(series: &Vsop87Series, t: f64) -> f64` - Series evaluator
  - ✅ `evaluate_vsop87_term(term: &Vsop87Term, t: f64) -> f64` - Individual term evaluator
  - ✅ `heliocentric_to_geocentric(...) -> Result<RaDec>` - Coordinate conversion
  - ✅ `get_planet_vsop87_data(planet: Planet) -> Option<PlanetVsop87Data>` - Data accessor
- ✅ Extend `CelestialObject` enum to include `Planet(Planet)` variant
- ✅ Define test cases
  - ✅ Test planet name parsing and enum functionality
  - ✅ Test VSOP87 term evaluation (cosine calculations)
  - ✅ Known planet positions at J2000.0 epoch (reference: JPL Horizons)
  - ✅ Accuracy requirements:
    - Inner planets (Mercury, Venus, Mars): ~1 arcminute (simplified) vs <1 arcsecond (full)
    - Outer planets (Jupiter, Saturn): ~1-2 arcminutes (simplified) vs <1 arcsecond (full)
    - Distant planets (Uranus, Neptune): ~2-5 arcminutes (simplified) vs <1 arcsecond (full)

**Implementation Details:**
- Created `src/planets.rs` module with complete data structure design
- VSOP87 series formula: `(L0 + L1×t + L2×t² + L3×t³ + L4×t⁴ + L5×t⁵) / 10^8`
- Each series term: `A × cos(B + C × t)` where t = Julian centuries from J2000.0
- Coordinate conversion: Heliocentric ecliptic (L, B, R) → Geocentric equatorial (RA, Dec)
- Integration with existing `CelestialObject` enum for unified API
- Placeholder implementations ready for Step B2 (data acquisition) and B3 (core implementation)

**Step B2: VSOP87 Data Acquisition** ✅ (Completed)
- ✅ Obtain VSOP87 coefficients
  - ✅ Researched official VSOP87 source (IMCCE FTP server)
  - ✅ Selected truncated/simplified version for initial implementation
  - ✅ Organized coefficients by planet and variable (L, B, R)
  - ✅ Implemented Mercury coefficients as proof of concept
  - ✅ Created placeholder functions for other planets
- ✅ Create coefficient storage structure
  - ✅ Designed efficient storage using const arrays (compile-time constants)
  - ✅ Implemented `Vsop87Term`, `Vsop87Series`, and `PlanetVsop87Data` structures
  - ✅ Created `get_planet_vsop87_data()` function with planet-specific data accessors
  - ✅ Documented coefficient format and source in code comments
- ✅ **📚 Educational Summary**: Added VSOP87 data structure explanation to OVERVIEW.md
  - ✅ Documented coefficient format and organization
  - ✅ Explained storage strategy (compile-time constants chosen for performance)
  - ✅ Added examples of coefficient data structure
  - ✅ Documented data acquisition process and sources

**Implementation Details:**
- Storage strategy: Compile-time constants (const arrays) for performance and simplicity
- Data structure: `Vsop87Term` (A, B, C) → `Vsop87Series` (series_0-5) → `PlanetVsop87Data` (L, B, R)
- Mercury implementation: Complete truncated coefficient set with ~10-20 terms per series
- Other planets: Placeholder structure ready for coefficient population
- Tests: 6 tests passing, verifying data structure and coefficient access
- Documentation: Comprehensive inline docs and OVERVIEW.md section on coefficient storage

**Step B3: Core Implementation - VSOP87 Evaluation** ⏳
- [ ] Implement VSOP87 series evaluator
  - Function to evaluate trigonometric series: `Σ(A * cos(B + C*t))`
  - Handle time argument (Julian centuries from J2000.0)
  - Optimize for performance (consider SIMD if beneficial)
  - Add logging for series evaluation (debug level)
- [ ] Implement L (longitude) calculation
  - Evaluate L0, L1, L2, L3, L4, L5 series
  - Combine series: `L = (L0 + L1*t + L2*t² + ...) / 10^8`
  - Convert to radians
- [ ] Implement B (latitude) calculation
  - Evaluate B0, B1, B2, B3, B4 series
  - Combine series similarly
- [ ] Implement R (radius) calculation
  - Evaluate R0, R1, R2, R3, R4 series
  - Combine series (in AU)
- [ ] **📚 Educational Summary**: Add VSOP87 evaluation algorithms to OVERVIEW.md
  - Document series evaluation algorithm
  - Explain time variable (Julian centuries) calculation
  - Add examples of series combination formula
  - Document performance considerations

**Step B4: Coordinate Conversion** ⏳
- [ ] Convert heliocentric ecliptic (L, B, R) to heliocentric equatorial
  - Apply ecliptic-to-equatorial rotation
  - Handle obliquity of the ecliptic (time-dependent)
- [ ] Convert to geocentric coordinates (optional, for apparent positions)
  - Account for Earth's position
  - Apply light-time correction (optional for initial version)
- [ ] Convert to RA/Dec
  - Final conversion to equatorial coordinates
  - Return `RaDec` structure
- [ ] **📚 Educational Summary**: Add coordinate conversion pipeline to OVERVIEW.md
  - Document heliocentric to geocentric conversion
  - Explain ecliptic to equatorial transformation
  - Add formulas for coordinate system conversions
  - Document light-time correction (if implemented)

**Step B5: Testing & Validation** ⏳
- [ ] Write unit tests for VSOP87 series evaluation
  - Test with known coefficient sets
  - Verify trigonometric calculations
  - Test edge cases (t = 0, large t values)
- [ ] Write unit tests for each planet
  - Test Mercury, Venus, Mars at J2000.0
  - Compare with reference values (JPL Horizons or Meeus)
  - Verify accuracy (within arcminutes for simplified, arcseconds for full)
- [ ] Write integration tests
  - Test multiple planets at same epoch
  - Test planets at different epochs (past and future)
  - Round-trip tests (if applicable)
- [ ] Performance testing
  - Benchmark planet position calculations
  - Ensure reasonable performance (< 10ms per planet)
- [ ] **📚 Educational Summary**: Add validation methodology to OVERVIEW.md
  - Document test strategy and reference sources
  - Explain accuracy requirements and validation approach
  - Add performance benchmarks and optimization notes

**Step B6: CLI Integration** ⏳
- [ ] Extend `position` command to support planets
  - Add planet names to object parsing
  - Update `CelestialObject` enum or create separate `Planet` enum
  - Integrate with existing position command structure
- [ ] Add new `planets` subcommand (optional)
  - List all available planets
  - Show planet positions for multiple planets
- [ ] Update help text and documentation
- [ ] Add example usage to readme
- [ ] **📚 Educational Summary**: Add CLI usage examples to OVERVIEW.md
  - Document planet position command usage
  - Add examples of planet calculations
  - Explain output format and interpretation

**Step B7: Logging & Error Handling** ⏳
- [ ] Add comprehensive logging
  - Debug: VSOP87 series values, intermediate calculations
  - Info: Planet calculations, coordinate conversions
  - Warn: Out-of-range time arguments, potential accuracy issues
- [ ] Enhance error handling
  - Validate time arguments (reasonable Julian Date range)
  - Handle missing coefficient data gracefully
  - Provide meaningful error messages for invalid planets
- [ ] **📚 Educational Summary**: Add error handling patterns to OVERVIEW.md
  - Document validation strategies
  - Explain error recovery approaches
  - Add examples of error handling for edge cases

---

#### Part C: Integration & Documentation

**Step C1: Module Organization** ⏳
- [ ] Review module structure
  - Consider if `planets.rs` should be separate from `celestial.rs`
  - Ensure proper module exports in `lib.rs`
  - Update project structure documentation
- [ ] Code organization and cleanup
  - Ensure consistent code style
  - Add comprehensive doc comments
  - Remove any placeholder code

**Step C2: Documentation Updates** ⏳
- [ ] Update readme.md
  - Add ECEF/ECI conversion examples
  - Add planet position examples
  - Update command documentation
  - Add references to VSOP87 and coordinate systems
- [ ] Add inline documentation
  - Document all public functions
  - Add examples to doc comments
  - Document coordinate system conventions
- [ ] **📚 Educational Summary**: Finalize OVERVIEW.md with complete VSOP87 documentation
  - Review and consolidate all VSOP87 sections
  - Add comprehensive examples and use cases
  - Document limitations and future enhancements
  - Create reference section for VSOP87 resources

**Step C3: Final Testing & Validation** ⏳
- [ ] Run full test suite
  - Ensure all existing tests still pass
  - Run new tests for ECEF/ECI and planets
  - Check test coverage
- [ ] Integration testing
  - Test complete workflows
  - Test CLI with various combinations
  - Verify output formatting
- [ ] Performance validation
  - Ensure no regressions in existing functionality
  - Verify new features perform acceptably

**Step C4: Changelog & Version Bump** ⏳
- [ ] Update version number (0.1 → 0.2)
- [ ] Document all changes in changelog
- [ ] Update feature list
- [ ] Tag release (if using version control)

---



## Part 2: Space Weather Web Service (Coming Next)

### Overview
A REST API service for fetching, caching, and serving space weather data critical for satellite operations. This service will complement the CLI tool by providing real-time and historical space weather information.

### Planned Features
- **Data Fetching**: Integration with NOAA Space Weather API and similar sources
- **Local Caching**: Reduce API calls and improve response times
- **REST Endpoints**: Clean API for querying current conditions and historical data
- **Data Storage**: Historical data storage in SQLite or PostgreSQL
- **Production Features**: Rate limiting and authentication
- **Deployment**: Self-hosted on personal server rack

### Use Cases
- Satellite operators monitoring space weather conditions
- Mission planning based on historical space weather patterns
- Real-time alerts for solar flares and geomagnetic storms
- Radiation level monitoring for space missions


### Testing Strategy

**Unit Tests:**
- Each function should have comprehensive unit tests
- Test edge cases, boundary conditions, and error conditions
- Use known reference values where possible
- Aim for >90% code coverage

**Integration Tests:**
- Test complete workflows end-to-end
- Test CLI commands with various inputs
- Verify output format and accuracy

**Validation Tests:**
- Compare results with authoritative sources (JPL Horizons, Meeus algorithms)
- Test round-trip conversions for accuracy
- Verify numerical stability

**Performance Tests:**
- Benchmark critical calculations
- Ensure no significant performance regressions
- Profile if performance issues arise

### Logging Strategy

**Log Levels:**
- **Error**: Calculation failures, invalid inputs
- **Warn**: Potential accuracy issues, out-of-range inputs
- **Info**: Major operations (coordinate transformations, planet calculations)
- **Debug**: Detailed intermediate calculations, matrix values, series evaluations

**Logging Points:**
- Function entry/exit for major operations
- Input validation results
- Intermediate calculation results (debug level)
- Transformation matrices (debug level)
- VSOP87 series evaluations (debug level)

### Still To Implement (Future Sessions - Part 1)
- More advanced orbital propagation
- Asteroid and comet positions
- Stellar positions and proper motion
- Advanced coordinate system transformations (GCRS, ITRS)
- VSOP87 planetary position calculations (Part B of Session 2)

### Part 2: Space Weather Web Service (Planned)
- REST API development
- NOAA Space Weather API integration
- Data caching and storage systems
- Rate limiting and authentication
- Server deployment and monitoring

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

Default location is Everett WA 47.90879667196036, -122.25035871889084

```bash
# Build
cargo build

# Test
cargo test

# Show all commands
cargo run -- --help
```

## Available Commands

### 1. `rise-set` - Calculate Rise/Set Times
Calculate when the Sun or Moon rises and sets at a specific location.

**Parameters:**
- `--object` or `-j`: Object name ("sun" or "moon")
- `--latitude` or `-a`: Observer latitude in degrees (positive = North, negative = South)
- `--longitude` or `-o`: Observer longitude in degrees (positive = East, negative = West)
- `--date` or `-d`: (Optional) Date in YYYY-MM-DD format (defaults to today)

**Examples:**
```bash
# Sunrise/sunset in New York today
cargo run -- rise-set --object "sun" --latitude 40.7128 --longitude=-74.0060

# Sunrise/sunset in Seattle on Christmas 2024
cargo run -- rise-set --object "sun" --latitude 47.6061 --longitude=-122.3328 --date "2024-12-25"

# Moonrise/moonset on summer solstice
cargo run -- rise-set --object "moon" --latitude 47.6061 --longitude=-122.3328 --date "2024-06-21"
```

### 2. `position` - Calculate Celestial Object Position
Calculate the Right Ascension and Declination of the Sun or Moon.

**Parameters:**
- `--object` or `-o`: Object name ("sun" or "moon")
- `--date` or `-d`: Date in YYYY-MM-DD format

**Examples:**
```bash
# Sun position on summer solstice
cargo run -- position --object "sun" --date "2024-06-21"
# Output: RA: 06:02:52, Dec: +23:26:02

# Moon position on a specific date
cargo run -- position --object "moon" --date "2024-01-01"
# Output: RA: 10:58:11, Dec: +10:02:27
```

### 3. `time` - Calculate Julian Date and Sidereal Time
Convert calendar dates to Julian Date and calculate Greenwich Mean Sidereal Time.

**Parameters:**
- `--date` or `-d`: Date in YYYY-MM-DD format
- `--time` or `-t`: (Optional) Time in HH:MM:SS format (defaults to 12:00:00)

**Examples:**
```bash
# Julian Date at noon
cargo run -- time --date "2024-01-01"
# Output: JD: 2460311.500000, GMST: 06:44:33

# Julian Date at specific time
cargo run -- time --date "2024-01-01" --time "18:30:45"
# Output: JD: 2460311.271354, GMST: 01:14:24
```

### 4. `convert` - Coordinate System Conversions
Convert between different coordinate systems: equatorial (RA/Dec), horizontal (Alt/Az), and Earth-centered (ECEF/ECI).

**Parameters:**
- `--from` or `-f`: Source coordinate system ("ra-dec", "alt-az", "ecef", or "eci")
- `--to` or `-t`: Target coordinate system ("alt-az", "ra-dec", "eci", or "ecef")
- `--coords` or `-c`: Coordinates to convert
- `--gmst`: (Optional) Greenwich Mean Sidereal Time in hours (0-24). If not provided, calculated from current time. Only used for ECEF/ECI conversions.

**Coordinate Formats:**
- RA/Dec: "hours,degrees" (e.g., "12.5,45.0")
- Alt/Az: "altitude,azimuth" in degrees (e.g., "45.0,180.0")
- ECEF/ECI: "x,y,z" in meters (e.g., "6378137.0,0.0,0.0")

**Examples:**
```bash
# Convert RA/Dec to Alt/Az (uses current time and default location: 47.9088°N, 122.2503°W)
cargo run -- convert --from ra-dec --to alt-az --coords "12.5,45.0"

# Convert Alt/Az to RA/Dec
cargo run -- convert --from alt-az --to ra-dec --coords "45.0,180.0"

# Convert ECEF to ECI (auto-calculates GMST from current time)
cargo run -- convert --from ecef --to eci --coords "6378137.0,0.0,0.0"

# Convert ECEF to ECI with specified GMST
cargo run -- convert --from ecef --to eci --coords "6378137.0,0.0,0.0" --gmst 12.5

# Convert ECI to ECEF (auto-calculates GMST from current time)
cargo run -- convert --from eci --to ecef --coords "6378137.0,0.0,0.0"
```

**Note:** 
- RA/Dec ↔ Alt/Az conversions vary based on current time and observer location (hardcoded as 47.9088°N, 122.2503°W).
- ECEF/ECI conversions require GMST (Greenwich Mean Sidereal Time). If `--gmst` is not provided, it is automatically calculated from the current time.
- ECEF/ECI coordinates are in meters. Typical values range from Earth's surface (~6,378,137 m) to geosynchronous orbit (~42,164,000 m).

### 5. `orbital` - Orbital Mechanics (Basic)
Display orbital elements (full state vector conversion coming soon).

**Parameters:**
- `--semi-major` or `-s`: Semi-major axis in km
- `--eccentricity` or `-e`: Orbital eccentricity (0-1)
- `--inclination` or `-i`: Orbital inclination in degrees

**Example:**
```bash
# ISS-like orbit
cargo run -- orbital --semi-major 6778 --eccentricity 0.0001 --inclination 51.6
# Output: Orbital: a=6778 km, e=0.0001, i=51.6°
```

## Formula Reference

### Time Systems
We utilize the full Julian day, vs the reduced. 
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

### ECEF/ECI Transformations (References)
- [ECEF](https://en.wikipedia.org/wiki/ECEF) - Earth-Centered Earth-Fixed coordinate system
- [ECI](https://en.wikipedia.org/wiki/Earth-centered_inertial) - Earth-Centered Inertial coordinate system
- [Sidereal Time](https://en.wikipedia.org/wiki/Sidereal_time) - Used for ECEF/ECI rotation

### Planetary Positions (References)
- [VSOP87](https://en.wikipedia.org/wiki/VSOP_(planets)) - Variable Speed Orbital Perturbations theory
- [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) - Ephemeris data for validation
- [Planetary Theory](https://en.wikipedia.org/wiki/Planetary_theory) - Mathematical models for planet positions

## Changelog

### Version 0.2 (Session 2 - In Progress)
**Completed Features:**
- ✅ ECEF ↔ ECI coordinate transformations (Part A complete)
- ✅ Extended CLI commands for ECEF/ECI conversions
- ✅ Enhanced logging and error handling
- ✅ Comprehensive test coverage for transformations

**Planned Features:**
- Planet position calculations using VSOP87 theory (Part B)
- Extended CLI commands for planet calculations

**Breaking Changes:**
- None expected

**Improvements:**
- Comprehensive test coverage for new features
- Performance optimizations for VSOP87 calculations
- Enhanced documentation and examples

---

### Version 0.3 (Planned - Part 2: Space Weather Service)
**Planned Features:**
- REST API for space weather data
- NOAA Space Weather API integration
- Data caching and historical storage
- Rate limiting and authentication
- Server deployment documentation

---

### Version 0.1 (Sessions 1-4 - Completed)
**Features:**
- Julian Date and Sidereal Time calculations
- Solar and lunar position calculations (RA/Dec)
- Rise/set times for Sun and Moon
- Coordinate conversions (RA/Dec ↔ Alt/Az)
- Kepler's equation solver
- Orbital period calculations
- Orbital elements to state vectors

**Phase 1 - Basic Structure:**
- Project structure with proper Rust organization
- CLI interface with subcommands
- Logging and error handling systems
- Module structure and documentation

**Phase 2 - Core Astronomy Calculations:**
- Julian Date calculation
- Greenwich Mean Sidereal Time (GMST)
- Local Sidereal Time (LST)
- Solar Position calculations
- Rise/Set Times calculations
- Comprehensive testing

**Phase 3 - Coordinate Conversions:**
- RA/Dec → Alt/Az conversions
- Alt/Az → RA/Dec conversions
- CLI Integration
- Real-time LST usage
- Comprehensive testing

**Phase 4 - Lunar & Orbital Mechanics:**
- Lunar Position calculations
- Lunar Rise/Set times
- Kepler's Equation solver
- Orbital Period calculations
- Elements → State Vectors conversion
- 29 comprehensive tests

