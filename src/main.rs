use clap::{Parser, Subcommand};
use cli_astro_calc::Result;
use log::{info, error};

/// A command-line astronomy and orbital mechanics calculator
#[derive(Parser)]
#[command(name = "cli-astro-calc")]
#[command(about = "A Swiss Army knife for space calculations")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    
    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert between coordinate systems (RA/Dec, Alt/Az, ECEF, ECI)
    Convert {
        /// Input coordinate system
        #[arg(short, long)]
        from: String,
        /// Output coordinate system  
        #[arg(short, long)]
        to: String,
        /// Input coordinates (format depends on system)
        #[arg(short, long)]
        coords: String,
    },
    /// Calculate rise/set times for celestial objects
    RiseSet {
        /// Celestial object name
        #[arg(short = 'j', long)]
        object: String,
        /// Observer latitude (degrees)
        #[arg(short = 'a', long)]
        latitude: f64,
        /// Observer longitude (degrees)
        #[arg(short = 'o', long)]
        longitude: f64,
    },
    /// Compute solar/lunar positions
    Position {
        /// Object: sun or moon
        #[arg(short, long)]
        object: String,
        /// Date (YYYY-MM-DD format)
        #[arg(short, long)]
        date: String,
    },
    /// Calculate Julian dates and sidereal time
    Time {
        /// Date (YYYY-MM-DD format)
        #[arg(short, long)]
        date: String,
        /// Time (HH:MM:SS format)
        #[arg(short, long)]
        time: Option<String>,
    },
    /// Convert orbital elements to position
    Orbital {
        /// Semi-major axis (km)
        #[arg(short, long)]
        semi_major: f64,
        /// Eccentricity
        #[arg(short, long)]
        eccentricity: f64,
        /// Inclination (degrees)
        #[arg(short, long)]
        inclination: f64,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    // Initialize logging
    init_logging(cli.verbose);
    
    info!("Starting CLI Astro Calc");
    
    match cli.command {
        Commands::Convert { from, to, coords } => {
            println!("Converting coordinates from {} to {}: {}", from, to, coords);
            
            // Parse the coordinate string and convert based on from/to types
            match (from.to_lowercase().as_str(), to.to_lowercase().as_str()) {
                ("ra-dec" | "radec", "alt-az" | "altaz") => {
                    // Parse RA/Dec format: "HH:MM:SS,DD:MM:SS" or "hours,degrees"
                    let result = parse_and_convert_radec_to_altaz(&coords)?;
                    println!("\nConversion Result:");
                    println!("Altitude:  {:.4}° ({:.2}° above horizon)", result.alt, result.alt);
                    println!("Azimuth:   {:.4}° (measured clockwise from North)", result.az);
                    
                    // Format as DMS
                    let alt_deg = result.alt.abs() as i32;
                    let alt_min = ((result.alt.abs() - alt_deg as f64) * 60.0) as i32;
                    let alt_sec = ((result.alt.abs() - alt_deg as f64 - alt_min as f64 / 60.0) * 3600.0) as i32;
                    let alt_sign = if result.alt >= 0.0 { "+" } else { "-" };
                    
                    let az_deg = result.az as i32;
                    let az_min = ((result.az - az_deg as f64) * 60.0) as i32;
                    let az_sec = ((result.az - az_deg as f64 - az_min as f64 / 60.0) * 3600.0) as i32;
                    
                    println!("\nFormatted:");
                    println!("Alt: {}{:02}°{:02}'{:02}\"", alt_sign, alt_deg, alt_min, alt_sec);
                    println!("Az:  {:03}°{:02}'{:02}\"", az_deg, az_min, az_sec);
                },
                ("alt-az" | "altaz", "ra-dec" | "radec") => {
                    // Parse Alt/Az format: "degrees,degrees"
                    let result = parse_and_convert_altaz_to_radec(&coords)?;
                    println!("\nConversion Result:");
                    println!("Right Ascension: {:.4} hours", result.ra);
                    println!("Declination:     {:.4}°", result.dec);
                    
                    // Format as HMS and DMS
                    let ra_h = result.ra as i32;
                    let ra_m = ((result.ra - ra_h as f64) * 60.0) as i32;
                    let ra_s = ((result.ra - ra_h as f64 - ra_m as f64 / 60.0) * 3600.0) as i32;
                    
                    let dec_deg = result.dec.abs() as i32;
                    let dec_min = ((result.dec.abs() - dec_deg as f64) * 60.0) as i32;
                    let dec_sec = ((result.dec.abs() - dec_deg as f64 - dec_min as f64 / 60.0) * 3600.0) as i32;
                    let dec_sign = if result.dec >= 0.0 { "+" } else { "-" };
                    
                    println!("\nFormatted:");
                    println!("RA:  {:02}:{:02}:{:02}", ra_h, ra_m, ra_s);
                    println!("Dec: {}{:02}°{:02}'{:02}\"", dec_sign, dec_deg, dec_min, dec_sec);
                },
                _ => {
                    return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                        format!("Unsupported conversion: {} to {}. Supported: ra-dec <-> alt-az", from, to)
                    ));
                }
            }
            
            info!("Coordinate conversion completed: {} -> {}", from, to);
        },
        Commands::RiseSet { object, latitude, longitude } => {
            println!("Calculating rise/set times for {} at ({}, {})", object, latitude, longitude);
            
            // Parse today's date (use current date if not specified)
            let date_time = chrono::Utc::now();
            
            // Create observer location
            let location = cli_astro_calc::celestial::ObserverLocation {
                latitude,
                longitude,
                elevation: 0.0, // Sea level by default
            };
            
            // Calculate rise/set times
            let rise_set = cli_astro_calc::celestial::calculate_rise_set_times(
                match object.to_lowercase().as_str() {
                    "sun" => cli_astro_calc::celestial::CelestialObject::Sun,
                    "moon" => cli_astro_calc::celestial::CelestialObject::Moon,
                    _ => return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                        format!("Unknown celestial object: {}", object)
                    )),
                },
                location,
                date_time
            )?;
            
            // Display results
            if let Some(rise) = rise_set.rise {
                println!("Rise time: {} UTC", rise.format("%Y-%m-%d %H:%M:%S"));
                println!("         : {} local", rise.format("%H:%M:%S"));
            } else {
                println!("Rise time: Object does not rise today (polar region)");
            }
            
            if let Some(set) = rise_set.set {
                println!("Set time:  {} UTC", set.format("%Y-%m-%d %H:%M:%S"));
                println!("         : {} local", set.format("%H:%M:%S"));
            } else {
                println!("Set time:  Object does not set today (polar region)");
            }
            
            info!("Rise/set calculation completed for {}", object);
        },
        Commands::Position { object, date } => {
            println!("Computing {} position for {}", object, date);
            
            // Parse the date
            let date_time = parse_date_time(&date, None)?;
            
            // Calculate position
            let position = cli_astro_calc::celestial::calculate_position(
                match object.to_lowercase().as_str() {
                    "sun" => cli_astro_calc::celestial::CelestialObject::Sun,
                    "moon" => cli_astro_calc::celestial::CelestialObject::Moon,
                    _ => return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                        format!("Unknown celestial object: {}", object)
                    )),
                },
                date_time
            )?;
            
            // Display results
            println!("Right Ascension: {:.4} hours", position.ra);
            println!("Declination: {:.4} degrees", position.dec);
            
            // Format RA as HH:MM:SS
            let ra_hours = position.ra as u32;
            let ra_minutes = ((position.ra - ra_hours as f64) * 60.0) as u32;
            let ra_seconds = ((position.ra - ra_hours as f64 - ra_minutes as f64 / 60.0) * 3600.0) as u32;
            println!("RA: {:02}:{:02}:{:02}", ra_hours, ra_minutes, ra_seconds);
            
            // Format Dec as DD:MM:SS
            let dec_degrees = position.dec.abs() as u32;
            let dec_minutes = ((position.dec.abs() - dec_degrees as f64) * 60.0) as u32;
            let dec_seconds = ((position.dec.abs() - dec_degrees as f64 - dec_minutes as f64 / 60.0) * 3600.0) as u32;
            let dec_sign = if position.dec >= 0.0 { "+" } else { "-" };
            println!("Dec: {}{:02}:{:02}:{:02}", dec_sign, dec_degrees, dec_minutes, dec_seconds);
            
            info!("Position calculation completed for {}", object);
        },
        Commands::Time { date, time } => {
            println!("Calculating Julian date and sidereal time for {}", date);
            
            // Parse the date and time
            let date_time = parse_date_time(&date, time.as_deref())?;
            
            // Calculate Julian Date
            let jd = cli_astro_calc::time::julian_date(date_time);
            println!("Julian Date: {:.6}", jd);
            
            // Calculate GMST
            let gmst = cli_astro_calc::time::greenwich_mean_sidereal_time(jd);
            println!("Greenwich Mean Sidereal Time: {:.4} hours", gmst);
            
            // Format GMST as HH:MM:SS
            let gmst_hours = gmst as u32;
            let gmst_minutes = ((gmst - gmst_hours as f64) * 60.0) as u32;
            let gmst_seconds = ((gmst - gmst_hours as f64 - gmst_minutes as f64 / 60.0) * 3600.0) as u32;
            println!("GMST: {:02}:{:02}:{:02}", gmst_hours, gmst_minutes, gmst_seconds);
            
            info!("Time calculation completed for {}", date);
        },
        Commands::Orbital { semi_major, eccentricity, inclination } => {
            println!("Converting orbital elements: a={}, e={}, i={}", semi_major, eccentricity, inclination);
            // TODO: Implement orbital conversion
            info!("Orbital conversion requested");
        },
    }
    
    info!("CLI Astro Calc completed successfully");
    Ok(())
}

fn init_logging(verbose: bool) {
    let log_level = if verbose { "debug" } else { "info" };
    std::env::set_var("RUST_LOG", log_level);
    env_logger::init();
}

/// Parse date and optional time string into DateTime<Utc>
fn parse_date_time(date_str: &str, time_str: Option<&str>) -> Result<chrono::DateTime<chrono::Utc>> {
    use chrono::{DateTime, Utc, NaiveDate, NaiveTime};
    
    // Parse date (YYYY-MM-DD format)
    let date = NaiveDate::parse_from_str(date_str, "%Y-%m-%d")
        .map_err(|e| cli_astro_calc::AstroError::InvalidTime(format!("Invalid date format: {}", e)))?;
    
    // Parse time if provided, otherwise default to noon
    let time = if let Some(time_str) = time_str {
        NaiveTime::parse_from_str(time_str, "%H:%M:%S")
            .map_err(|e| cli_astro_calc::AstroError::InvalidTime(format!("Invalid time format: {}", e)))?
    } else {
        NaiveTime::from_hms_opt(12, 0, 0).unwrap() // Default to noon
    };
    
    // Combine date and time, then convert to UTC
    let naive_datetime = date.and_time(time);
    let datetime = DateTime::from_naive_utc_and_offset(naive_datetime, Utc);
    
    Ok(datetime)
}

/// Parse RA/Dec coordinates and convert to Alt/Az
/// Format: "hours,degrees" or "HH:MM:SS,DD:MM:SS"
/// Requires observer location and current time (uses current time and assumes New York for demo)
fn parse_and_convert_radec_to_altaz(coords: &str) -> Result<cli_astro_calc::coordinates::AltAz> {
    use cli_astro_calc::coordinates::{RaDec, ra_dec_to_alt_az};
    use cli_astro_calc::time::{julian_date, greenwich_mean_sidereal_time, local_sidereal_time};
    
    // Parse coordinates (simple format: "hours,degrees")
    let parts: Vec<&str> = coords.split(',').collect();
    if parts.len() != 2 {
        return Err(cli_astro_calc::AstroError::InvalidCoordinate(
            "Expected format: hours,degrees (e.g., '12.5,45.0')".to_string()
        ));
    }
    
    let ra: f64 = parts[0].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid RA value".to_string()))?;
    let dec: f64 = parts[1].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid Dec value".to_string()))?;
    
    let ra_dec = RaDec { ra, dec };
    
    // Use current time and New York coordinates for demonstration
    let now = chrono::Utc::now();
    let jd = julian_date(now);
    let gmst = greenwich_mean_sidereal_time(jd);
    let observer_lat = 40.7128; // New York
    let observer_lon = -74.0060;
    let lst = local_sidereal_time(gmst, observer_lon);
    
    println!("Using current LST: {:.4} hours at New York (lat={}, lon={})", lst, observer_lat, observer_lon);
    
    ra_dec_to_alt_az(ra_dec, observer_lat, observer_lon, lst)
}

/// Parse Alt/Az coordinates and convert to RA/Dec
/// Format: "altitude,azimuth" in degrees
/// Requires observer location and current time (uses current time and assumes New York for demo)
fn parse_and_convert_altaz_to_radec(coords: &str) -> Result<cli_astro_calc::coordinates::RaDec> {
    use cli_astro_calc::coordinates::{AltAz, alt_az_to_ra_dec};
    use cli_astro_calc::time::{julian_date, greenwich_mean_sidereal_time, local_sidereal_time};
    
    // Parse coordinates (format: "altitude,azimuth")
    let parts: Vec<&str> = coords.split(',').collect();
    if parts.len() != 2 {
        return Err(cli_astro_calc::AstroError::InvalidCoordinate(
            "Expected format: altitude,azimuth in degrees (e.g., '45.0,180.0')".to_string()
        ));
    }
    
    let alt: f64 = parts[0].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid altitude value".to_string()))?;
    let az: f64 = parts[1].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid azimuth value".to_string()))?;
    
    let alt_az = AltAz { alt, az };
    
    // Use current time and New York coordinates for demonstration
    let now = chrono::Utc::now();
    let jd = julian_date(now);
    let gmst = greenwich_mean_sidereal_time(jd);
    let observer_lat = 40.7128; // New York
    let observer_lon = -74.0060;
    let lst = local_sidereal_time(gmst, observer_lon);
    
    println!("Using current LST: {:.4} hours at New York (lat={}, lon={})", lst, observer_lat, observer_lon);
    
    alt_az_to_ra_dec(alt_az, observer_lat, observer_lon, lst)
}
