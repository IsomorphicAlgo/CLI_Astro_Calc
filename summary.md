# CLI Astro Calc - Project Architecture Summary

## Overview
This Rust project is a command-line astronomy calculator that demonstrates Rust's memory safety, ownership system, and modular architecture. It's designed to teach both Rust programming concepts and basic astrophysics calculations.

## File Structure and Interactions

### 1. **Cargo.toml** - Project Configuration
- **Purpose**: Defines project metadata, dependencies, and build configuration
- **Rust Concepts**: Package management, dependency resolution
- **Key Dependencies**:
  - `clap`: Command-line argument parsing (uses derive macros)
  - `chrono`: Date/time handling (demonstrates external crate usage)
  - `log` + `env_logger`: Logging system (shows trait-based design)
  - `anyhow` + `thiserror`: Error handling (demonstrates custom error types)

### 2. **src/main.rs** - Application Entry Point
- **Purpose**: CLI interface and command routing
- **Rust Concepts Demonstrated**:
  - **Ownership**: CLI arguments are moved into the `Cli` struct
  - **Pattern Matching**: `match` expressions for command handling
  - **Error Propagation**: `Result<()>` return type with `?` operator
  - **Lifetimes**: Implicit lifetime management in string handling

```rust
// Ownership example: clap takes ownership of command line args
let cli = Cli::parse(); // Cli owns the parsed arguments

// Pattern matching with ownership transfer
match cli.command {
    Commands::RiseSet { object, latitude, longitude } => {
        // object, latitude, longitude are moved here
        // They're now owned by this scope
    }
}
```

### 3. **src/lib.rs** - Library Root and Module Organization
- **Purpose**: Defines the public API and module structure
- **Rust Concepts**:
  - **Module System**: `pub mod` declarations expose modules
  - **Re-exports**: `pub use` makes types available at crate root
  - **Error Types**: Custom error enum with `thiserror` derive macro

```rust
// Module organization - each module is a separate file
pub mod coordinates;  // Coordinate system conversions
pub mod celestial;    // Celestial object calculations  
pub mod time;         // Time calculations
pub mod orbital;      // Orbital mechanics

// Re-export for convenience
pub use error::{AstroError, Result};
```

### 4. **src/coordinates.rs** - Coordinate System Conversions
- **Purpose**: Convert between different astronomical coordinate systems
- **Rust Concepts**:
  - **Struct Definition**: `#[derive(Debug, Clone, Copy)]` for data types
  - **Copy Semantics**: All coordinate structs implement `Copy` (stack-allocated)
  - **Ownership**: Functions take ownership of coordinate structs (but Copy means no actual move)

```rust
// Copy semantics - these structs are copied, not moved
#[derive(Debug, Clone, Copy)]
pub struct RaDec {
    pub ra: f64,   // Right Ascension in hours
    pub dec: f64,  // Declination in degrees
}

// Function signature shows ownership transfer (but Copy prevents actual move)
pub fn ra_dec_to_alt_az(ra_dec: RaDec, observer_lat: f64, observer_lon: f64, lst: f64) -> Result<AltAz>
```

**Astrophysics Context**: 
- **RA/Dec**: Equatorial coordinates (like latitude/longitude on Earth, but for the sky)
- **Alt/Az**: Horizontal coordinates (altitude/azimuth from observer's perspective)
- **ECEF/ECI**: Earth-centered coordinate systems for satellite calculations

### 5. **src/celestial.rs** - Celestial Object Calculations
- **Purpose**: Calculate positions and rise/set times for Sun and Moon
- **Rust Concepts**:
  - **Enums**: `CelestialObject` enum for type safety
  - **Option Types**: `RiseSetTimes` uses `Option<DateTime<Utc>>` for cases where objects don't rise/set
  - **External Dependencies**: Uses `chrono::DateTime<Utc>` (owned type)

```rust
// Enum for type safety - prevents invalid celestial objects
#[derive(Debug, Clone, Copy)]
pub enum CelestialObject {
    Sun,
    Moon,
}

// Option type handles cases where celestial objects don't rise/set (polar regions)
#[derive(Debug, Clone, Copy)]
pub struct RiseSetTimes {
    pub rise: Option<DateTime<Utc>>,  // Option owns the DateTime
    pub set: Option<DateTime<Utc>>,    // None if object doesn't rise/set
}
```

**Astrophysics Context**:
- **Rise/Set Times**: When celestial objects appear/disappear from observer's horizon
- **Observer Location**: Latitude/longitude affect what's visible and when
- **Polar Regions**: Some objects never rise or set depending on observer location

### 6. **src/time.rs** - Time Calculations
- **Purpose**: Julian dates and sidereal time calculations
- **Rust Concepts**:
  - **Pure Functions**: All functions are pure (no side effects, same input = same output)
  - **Ownership**: Functions take ownership of `DateTime<Utc>` but return owned `f64`
  - **No Mutable State**: All calculations are immutable

```rust
// Pure function - takes ownership, returns owned value
pub fn julian_date(date_time: DateTime<Utc>) -> f64 {
    // DateTime<Utc> is moved into this function
    // Returns owned f64 (Copy type)
}
```

**Astrophysics Context**:
- **Julian Date**: Continuous day count since January 1, 4713 BC (astronomical standard)
- **Sidereal Time**: "Star time" - Earth's rotation relative to stars, not Sun
- **GMST/LST**: Greenwich Mean Sidereal Time and Local Sidereal Time

### 7. **src/orbital.rs** - Orbital Mechanics
- **Purpose**: Convert orbital elements to position/velocity vectors
- **Rust Concepts**:
  - **Array Types**: `[f64; 3]` for 3D vectors (stack-allocated, Copy)
  - **Struct Composition**: Complex structs containing arrays
  - **Ownership**: Functions take ownership of `OrbitalElements` struct

```rust
// Array type - stack allocated, implements Copy
#[derive(Debug, Clone, Copy)]
pub struct StateVector {
    pub position: [f64; 3],  // 3D position vector
    pub velocity: [f64; 3],  // 3D velocity vector
}
```

**Astrophysics Context**:
- **Orbital Elements**: Six parameters that completely describe an orbit
- **State Vectors**: Position and velocity at a specific time
- **Kepler's Laws**: Mathematical relationships governing orbital motion

## Memory Ownership Patterns

### 1. **Copy Types** (Stack Allocation)
```rust
// These types implement Copy - they're copied, not moved
f64, i32, bool, char
RaDec, AltAz, Ecef, Eci  // Our coordinate structs
CelestialObject          // Our enum
```

### 2. **Owned Types** (Heap Allocation)
```rust
// These types are moved, not copied
String                   // Dynamic string
Vec<T>                   // Dynamic vector
DateTime<Utc>           // From chrono crate
Box<T>                  // Heap-allocated box
```

### 3. **Borrowing Patterns**
```rust
// Functions can borrow instead of taking ownership
fn calculate_position(object: &CelestialObject, date: &DateTime<Utc>) -> Result<RaDec>
//                                 ^borrow        ^borrow
```

## Error Handling Strategy

### Custom Error Types
```rust
#[derive(Error, Debug)]
pub enum AstroError {
    #[error("Invalid coordinate value: {0}")]
    InvalidCoordinate(String),  // Owned String for error message
    
    #[error("Calculation error: {0}")]
    CalculationError(String),
    
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),  // Automatic conversion from std::io::Error
}
```

### Result Propagation
```rust
// Functions return Result<T, AstroError>
pub fn ra_dec_to_alt_az(...) -> Result<AltAz> {
    // Result<AltAz, AstroError> - AltAz is Copy, AstroError is owned
}

// Main function uses ? operator for error propagation
fn main() -> Result<()> {
    let cli = Cli::parse();  // Could fail, but clap handles it
    // ... other operations that could fail
    Ok(())  // Success case
}
```

## Module Interaction Flow

1. **main.rs** parses CLI arguments and owns them
2. **main.rs** calls functions from other modules, passing ownership
3. **lib.rs** re-exports types for easy access
4. **Individual modules** contain domain-specific logic
5. **Error types** flow back up through Result types
6. **All calculations** are pure functions with clear ownership

## Learning Points for Rust Beginners

### 1. **Ownership is Everywhere**
- Every value has exactly one owner
- Ownership can be transferred (moved) or borrowed
- Copy types avoid ownership issues for simple data

### 2. **Memory Safety Without Garbage Collection**
- Rust prevents use-after-free, double-free, and memory leaks at compile time
- No runtime overhead for memory management
- Stack allocation is preferred (Copy types)

### 3. **Error Handling is Explicit**
- No exceptions - all errors are handled through Result types
- `?` operator makes error propagation concise
- Custom error types provide meaningful error messages

### 4. **Module System Encourages Organization**
- Each file is a module with a specific purpose
- `pub` keyword controls what's exposed
- Dependencies flow in one direction (main.rs → lib.rs → modules)

## Astrophysics Learning Points

### 1. **Coordinate Systems**
- Different coordinate systems for different purposes
- Conversions require mathematical transformations
- Observer location affects coordinate values

### 2. **Time Systems**
- Astronomical time is different from civil time
- Julian dates provide continuous time reference
- Sidereal time accounts for Earth's rotation

### 3. **Celestial Mechanics**
- Objects in space follow predictable patterns
- Mathematical models can predict positions
- Observer location affects what's visible

This architecture demonstrates how Rust's ownership system enables safe, efficient code while teaching fundamental astronomical concepts through a practical, working application.
