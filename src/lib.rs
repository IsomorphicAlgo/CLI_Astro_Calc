//! CLI Astro Calc - A command-line astronomy and orbital mechanics calculator
//! 
//! This library provides functions for various astronomical calculations including:
//! - Coordinate system conversions (RA/Dec, Alt/Az, ECEF, ECI)
//! - Rise/set times for celestial objects
//! - Solar/lunar position calculations
//! - Julian dates and sidereal time
//! - Basic orbital elements to position conversions

pub mod coordinates;
pub mod celestial;
pub mod time;
pub mod orbital;

/// Error types for the astronomy calculator
pub mod error {
    use thiserror::Error;

    #[derive(Error, Debug)]
    pub enum AstroError {
        #[error("Invalid coordinate value: {0}")]
        InvalidCoordinate(String),
        
        #[error("Invalid time value: {0}")]
        InvalidTime(String),
        
        #[error("Calculation error: {0}")]
        CalculationError(String),
        
        #[error("IO error: {0}")]
        IoError(#[from] std::io::Error),
    }

    pub type Result<T> = std::result::Result<T, AstroError>;
}

pub use error::{AstroError, Result};
