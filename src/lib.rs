pub mod coordinates;
pub mod celestial;
pub mod time;
pub mod orbital;

pub mod error {
    use thiserror::Error;

    #[derive(Error, Debug)]
    pub enum AstroError {
        #[error("Invalid coordinate: {0}")]
        InvalidCoordinate(String),
        
        #[error("Invalid time: {0}")]
        InvalidTime(String),
        
        #[error("Calculation error: {0}")]
        CalculationError(String),
        
        #[error("IO error: {0}")]
        IoError(#[from] std::io::Error),
    }

    pub type Result<T> = std::result::Result<T, AstroError>;
}

pub use error::{AstroError, Result};
