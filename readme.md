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
Convert between equatorial (RA/Dec) and horizontal (Alt/Az) coordinate systems.

**Parameters:**
- `--from` or `-f`: Source coordinate system ("ra-dec" or "alt-az")
- `--to` or `-t`: Target coordinate system ("alt-az" or "ra-dec")
- `--coords` or `-c`: Coordinates to convert

**Coordinate Formats:**
- RA/Dec: "hours,degrees" (e.g., "12.5,45.0")
- Alt/Az: "altitude,azimuth" in degrees (e.g., "45.0,180.0")

**Examples:**
```bash
# Convert RA/Dec to Alt/Az (uses current time and New York location)
cargo run -- convert --from ra-dec --to alt-az --coords "12.5,45.0"
# Output: Alt: +32°36'06", Az: 056°09'08"

# Convert Alt/Az to RA/Dec
cargo run -- convert --from alt-az --to ra-dec --coords "45.0,180.0"
# Output: RA: 07:03:25, Dec: -04°17'13"
```

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

