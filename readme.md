# CLI Astro Calc

A command-line astronomy and orbital mechanics calculator built in Rust. The purpose of this is to learn the Rust language whilst creating a custom tool for myself.

## Project Status: Phase 3 - Coordinate Conversions ✅

This project is being developed iteratively, with each phase building upon the previous one.

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

### Still To Implement
- Lunar position calculations
- Basic orbital elements to position conversions
- ECEF ↔ ECI transformations

## Project Structure

```
CLI_Astro_Calc/
├── Cargo.toml          # Project configuration and dependencies
├── src/
│   ├── main.rs         # CLI entry point and command parsing
│   ├── lib.rs          # Library structure and error types
│   ├── coordinates.rs  # Coordinate system conversions
│   ├── celestial.rs    # Celestial object calculations
│   ├── time.rs         # Time calculations (Julian dates, sidereal time)
│   └── orbital.rs      # Orbital mechanics calculations
├── tests/              # Integration tests
├── benches/            # Performance benchmarks
├── docs/               # Documentation
└── readme.md           # This file
```

## Dependencies

- **clap**: Command-line argument parsing
- **log**: Logging framework
- **env_logger**: Logging implementation
- **anyhow**: Error handling
- **thiserror**: Custom error types
- **chrono**: Date and time handling
- **serde**: Serialization support
- **criterion**: Benchmarking (dev dependency)

## Building and Running

### Prerequisites
- Rust 1.70+ (edition 2021)
- Cargo package manager

### Build
```bash
cargo build
```

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
# Output: Julian Date: 2460311.500000, GMST: 06:44:33

# Convert RA/Dec to Alt/Az (uses current time and New York location)
cargo run -- convert --from ra-dec --to alt-az --coords "12.5,45.0"
# Output: Alt: +32°36'06", Az: 056°09'08"

# Convert Alt/Az to RA/Dec (uses current time and New York location)
cargo run -- convert --from alt-az --to ra-dec --coords "45.0,180.0"
# Output: RA: 07:03:25, Dec: -04°17'13"

# Placeholder examples (not yet implemented):
cargo run -- position --object "moon" --date "2024-01-01"
cargo run -- orbital --semi-major 7000 --eccentricity 0.1 --inclination 45

# Enable verbose logging
cargo run -- --verbose rise-set --object "sun" --latitude 40.7128 --longitude=-74.0060

# Note: For negative values, use the = syntax: --longitude=-74.0060
```

### Test
```bash
cargo test
```

### Benchmark
```bash
cargo bench
```

## Function Descriptions

### Coordinate Conversions
- **RA/Dec to Alt/Az**: Convert from equatorial coordinates (Right Ascension/Declination) to horizontal coordinates (Altitude/Azimuth) for a given observer location and time
- **Alt/Az to RA/Dec**: Convert from horizontal coordinates to equatorial coordinates
- **ECEF/ECI**: Convert between Earth-Centered Earth-Fixed and Earth-Centered Inertial coordinate systems

### Celestial Calculations
- **Rise/Set Times**: Calculate when celestial objects rise and set for a given observer location
- **Solar/Lunar Positions**: Compute the current position of the Sun or Moon in the sky

### Time Calculations
- **Julian Dates**: Convert calendar dates to Julian Date system used in astronomy
- **Sidereal Time**: Calculate sidereal time (star time) for astronomical observations

### Orbital Mechanics
- **Orbital Elements to Position**: Convert classical orbital elements to position and velocity vectors
- **Orbital Period**: Calculate the period of an orbit from its semi-major axis

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

## Next Phase
Future phases will implement:
- Lunar position calculations (more complex than solar due to orbital perturbations)
- Orbital mechanics (elements to state vectors, Kepler's equations)
- ECEF ↔ ECI transformations (Earth-fixed to inertial frames)
- More celestial objects (planets, stars from catalogs)

---

**Important**: Remember to commit and push to GitHub after each phase completion!

