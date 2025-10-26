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

# Calculate sunrise/sunset times for today
cargo run -- rise-set --object "sun" --latitude 40.7128 --longitude=-74.0060

# For a specific date (Seattle on Christmas)
cargo run -- rise-set --object "sun" --latitude 47.6061 --longitude=-122.3328 --date "2024-12-25"

# Moon rise/set (note: west longitudes are negative)
cargo run -- rise-set --object "moon" --latitude 47.6061 --longitude=-122.3328 --date "2024-06-21"

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

