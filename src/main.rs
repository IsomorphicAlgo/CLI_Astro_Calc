use clap::{Parser, Subcommand};
use cli_astro_calc::Result;

/// A command-line astronomy and orbital mechanics calculator
#[derive(Parser)]
#[command(name = "cli-astro-calc")]
#[command(about = "An Astrological Calculator and exploration into RUST")]
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
    Convert {
        #[arg(short, long)]
        from: String,
        #[arg(short, long)]
        to: String,
        #[arg(short, long)]
        coords: String,
    },
    RiseSet {
        #[arg(short = 'j', long)]
        object: String,
        #[arg(short = 'a', long)]
        latitude: f64,
        #[arg(short = 'o', long)]
        longitude: f64,
        #[arg(short, long)]
        date: Option<String>,
    },
    Position {
        #[arg(short, long)]
        object: String,
        #[arg(short, long)]
        date: String,
    },
    Time {
        #[arg(short, long)]
        date: String,
        #[arg(short, long)]
        time: Option<String>,
    },
    Orbital {
        #[arg(short, long)]
        semi_major: f64,
        #[arg(short, long)]
        eccentricity: f64,
        #[arg(short, long)]
        inclination: f64,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.verbose);
    
    match cli.command {
        Commands::Convert { from, to, coords } => {
            match (from.to_lowercase().as_str(), to.to_lowercase().as_str()) {
                ("ra-dec" | "radec", "alt-az" | "altaz") => {
                    let result = parse_and_convert_radec_to_altaz(&coords)?;
                    let (alt_deg, alt_min, alt_sec, alt_sign) = format_angle(result.alt);
                    let (az_deg, az_min, az_sec, _) = format_angle(result.az);
                    println!("Alt: {}{:02}°{:02}'{:02}\"", alt_sign, alt_deg, alt_min, alt_sec);
                    println!("Az:  {:03}°{:02}'{:02}\"", az_deg, az_min, az_sec);
                },
                ("alt-az" | "altaz", "ra-dec" | "radec") => {
                    let result = parse_and_convert_altaz_to_radec(&coords)?;
                    let (ra_h, ra_m, ra_s) = format_time(result.ra);
                    let (dec_deg, dec_min, dec_sec, dec_sign) = format_angle(result.dec);
                    println!("RA:  {:02}:{:02}:{:02}", ra_h, ra_m, ra_s);
                    println!("Dec: {}{:02}°{:02}'{:02}\"", dec_sign, dec_deg, dec_min, dec_sec);
                },
                _ => {
                    return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                        format!("Unsupported conversion: {} to {}", from, to)
                    ));
                }
            }
        },
        Commands::RiseSet { object, latitude, longitude, date } => {
            let date_time = if let Some(date_str) = date {
                parse_date_time(&date_str, None)?
            } else {
                chrono::Utc::now()
            };
            
            let location = cli_astro_calc::celestial::ObserverLocation {
                latitude,
                longitude,
                elevation: 0.0,
            };
            
            let obj = match object.to_lowercase().as_str() {
                "sun" => cli_astro_calc::celestial::CelestialObject::Sun,
                "moon" => cli_astro_calc::celestial::CelestialObject::Moon,
                _ => return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                    format!("Unknown object: {}", object)
                )),
            };
            
            let rise_set = cli_astro_calc::celestial::calculate_rise_set_times(obj, location, date_time)?;
            
            match rise_set.rise {
                Some(t) => println!("Rise: {}", t.format("%H:%M:%S UTC")),
                None => println!("Rise: Does not rise"),
            }
            match rise_set.set {
                Some(t) => println!("Set:  {}", t.format("%H:%M:%S UTC")),
                None => println!("Set:  Does not set"),
            }
        },
        Commands::Position { object, date } => {
            let date_time = parse_date_time(&date, None)?;
            let obj = match object.to_lowercase().as_str() {
                "sun" => cli_astro_calc::celestial::CelestialObject::Sun,
                "moon" => cli_astro_calc::celestial::CelestialObject::Moon,
                _ => return Err(cli_astro_calc::AstroError::InvalidCoordinate(
                    format!("Unknown object: {}", object)
                )),
            };
            
            let pos = cli_astro_calc::celestial::calculate_position(obj, date_time)?;
            let (ra_h, ra_m, ra_s) = format_time(pos.ra);
            let (dec_deg, dec_min, dec_sec, dec_sign) = format_angle(pos.dec);
            
            println!("RA:  {:02}:{:02}:{:02}", ra_h, ra_m, ra_s);
            println!("Dec: {}{:02}°{:02}'{:02}\"", dec_sign, dec_deg, dec_min, dec_sec);
        },
        Commands::Time { date, time } => {
            let date_time = parse_date_time(&date, time.as_deref())?;
            let jd = cli_astro_calc::time::julian_date(date_time);
            let gmst = cli_astro_calc::time::greenwich_mean_sidereal_time(jd);
            let (h, m, s) = format_time(gmst);
            
            println!("JD:   {:.6}", jd);
            println!("GMST: {:02}:{:02}:{:02}", h, m, s);
        },
        Commands::Orbital { semi_major, eccentricity, inclination } => {
            println!("Orbital: a={} km, e={}, i={}°", semi_major, eccentricity, inclination);
        },
    }
    
    Ok(())
}

fn init_logging(verbose: bool) {
    let log_level = if verbose { "debug" } else { "info" };
    std::env::set_var("RUST_LOG", log_level);
    env_logger::init();
}

fn parse_date_time(date_str: &str, time_str: Option<&str>) -> Result<chrono::DateTime<chrono::Utc>> {
    use chrono::{DateTime, Utc, NaiveDate, NaiveTime};
    
    let date = NaiveDate::parse_from_str(date_str, "%Y-%m-%d")
        .map_err(|e| cli_astro_calc::AstroError::InvalidTime(format!("Invalid date: {}", e)))?;
    
    let time = if let Some(ts) = time_str {
        NaiveTime::parse_from_str(ts, "%H:%M:%S")
            .map_err(|e| cli_astro_calc::AstroError::InvalidTime(format!("Invalid time: {}", e)))?
    } else {
        NaiveTime::from_hms_opt(12, 0, 0).unwrap()
    };
    
    Ok(DateTime::from_naive_utc_and_offset(date.and_time(time), Utc))
}

fn format_time(hours: f64) -> (i32, i32, i32) {
    let h = hours as i32;
    let m = ((hours - h as f64) * 60.0) as i32;
    let s = ((hours - h as f64 - m as f64 / 60.0) * 3600.0) as i32;
    (h, m, s)
}

fn format_angle(degrees: f64) -> (i32, i32, i32, &'static str) {
    let deg = degrees.abs() as i32;
    let min = ((degrees.abs() - deg as f64) * 60.0) as i32;
    let sec = ((degrees.abs() - deg as f64 - min as f64 / 60.0) * 3600.0) as i32;
    let sign = if degrees >= 0.0 { "+" } else { "-" };
    (deg, min, sec, sign)
}

fn parse_and_convert_radec_to_altaz(coords: &str) -> Result<cli_astro_calc::coordinates::AltAz> {
    use cli_astro_calc::coordinates::{RaDec, ra_dec_to_alt_az};
    use cli_astro_calc::time::{julian_date, greenwich_mean_sidereal_time, local_sidereal_time};
    
    let parts: Vec<&str> = coords.split(',').collect();
    if parts.len() != 2 {
        return Err(cli_astro_calc::AstroError::InvalidCoordinate("Expected: hours,degrees".to_string()));
    }
    
    let ra: f64 = parts[0].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid RA".to_string()))?;
    let dec: f64 = parts[1].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid Dec".to_string()))?;
    
    let jd = julian_date(chrono::Utc::now());
    let gmst = greenwich_mean_sidereal_time(jd);
    let (lat, lon) = (40.7128, -74.0060);
    let lst = local_sidereal_time(gmst, lon);
    
    ra_dec_to_alt_az(RaDec { ra, dec }, lat, lon, lst)
}

fn parse_and_convert_altaz_to_radec(coords: &str) -> Result<cli_astro_calc::coordinates::RaDec> {
    use cli_astro_calc::coordinates::{AltAz, alt_az_to_ra_dec};
    use cli_astro_calc::time::{julian_date, greenwich_mean_sidereal_time, local_sidereal_time};
    
    let parts: Vec<&str> = coords.split(',').collect();
    if parts.len() != 2 {
        return Err(cli_astro_calc::AstroError::InvalidCoordinate("Expected: altitude,azimuth".to_string()));
    }
    
    let alt: f64 = parts[0].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid altitude".to_string()))?;
    let az: f64 = parts[1].trim().parse()
        .map_err(|_| cli_astro_calc::AstroError::InvalidCoordinate("Invalid azimuth".to_string()))?;
    
    let jd = julian_date(chrono::Utc::now());
    let gmst = greenwich_mean_sidereal_time(jd);
    let (lat, lon) = (40.7128, -74.0060);
    let lst = local_sidereal_time(gmst, lon);
    
    alt_az_to_ra_dec(AltAz { alt, az }, lat, lon, lst)
}
