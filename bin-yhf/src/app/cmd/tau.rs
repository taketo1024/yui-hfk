use log::error;
use yui_hfk::invariants::tau::tau_with_config;
use crate::app::{AppArgs, ConfigArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// Use the mirror of the knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Expected genus; aborts if it disagrees with the computed value
    #[arg(short, long)]
    pub g: Option<usize>,

    /// Filter candidates using the 4-genus bound
    #[arg(short='G', long)]
    pub g4: Option<usize>,

    #[command(flatten)]
    pub config: ConfigArgs,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 { self.log }
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    validate(args);

    let gd = load_diagram(&args.knot)?;
    let gd = if args.mirror { gd.rotate90() } else { gd };
    gd.log_diagram();

    let t = tau_with_config(&gd, args.g, args.g4, args.config.config());

    match t { 
        Ok(t) => Ok(format!("{t}")),
        Err(e) => {
            error!("Unable to compute τ.");
            Err(format!("{e:?}").into())
        }
    }
}

fn validate(args: &Args) { 
    if let Some(g) = args.g && let Some(g4) = args.g4 { 
        assert!(g4 <= g, "Invalid inequality: g4 = {g4} > g = {g}");
    }
}