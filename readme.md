# CLI Astro Calc

A command-line astronomy and orbital mechanics calculator built in Rust. The purpose of this is to learn the Rust language whilst creating a custom tool for myself.

## Project Status: Phase 1 - Basic Structure ✅

This project is being developed iteratively, with each phase building upon the previous one.

### Current Features (Phase 1)
- ✅ Basic project structure with proper Rust organization
- ✅ CLI interface using clap with subcommands for all planned features
- ✅ Logging system with verbose output option
- ✅ Error handling with custom error types
- ✅ Module structure for different calculation types
- ✅ Comprehensive documentation and type definitions

### Planned Features
- Convert between coordinate systems (RA/Dec, Alt/Az, ECEF, ECI)
- Calculate rise/set times for celestial objects from a given location
- Compute solar/lunar positions
- Calculate Julian dates and sidereal time
- Basic orbital elements to position conversions

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

# Example commands (currently show placeholder output)
cargo run -- convert --from "ra-dec" --to "alt-az" --coords "12:34:56,45:30:15"
cargo run -- rise-set --object "sun" --latitude 40.7128 --longitude=-74.0060
cargo run -- position --object "moon" --date "2024-01-01"
cargo run -- time --date "2024-01-01" --time "12:00:00"
cargo run -- orbital --semi-major 7000 --eccentricity 0.1 --inclination 45

# Enable verbose logging
cargo run -- --verbose convert --from "ra-dec" --to "alt-az" --coords "12:34:56,45:30:15"

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

## Next Phase
The next phase will implement the actual calculation algorithms, starting with the most fundamental functions like Julian Date calculations and basic coordinate conversions.

---

**Important**: Remember to commit and push to GitHub after each phase completion!

