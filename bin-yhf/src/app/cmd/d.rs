use log::error;
use yui_core::num::Ratio;
use yui_hfk::G0CalcConfig;
use yui_hfk::invariants::d::d_with_config;
use crate::app::{AppArgs, ConfigArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// The surgery coefficient p/q ≠ 0, e.g. "5", "5/2" or "-3/2"
    #[arg(allow_hyphen_values = true)]
    pub slope: Ratio<i64>,

    /// Use the mirror of the knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Expected genus; aborts if it disagrees with the computed value
    #[arg(short, long)]
    pub g: Option<usize>,

    /// Filter candidates using the 4-genus bound
    #[arg(short='G', long)]
    pub g4: Option<usize>,

    /// Filter candidates using a known value of the tau invariant
    #[arg(short, long)]
    pub tau: Option<isize>,

    #[command(flatten)]
    pub config: ConfigArgs,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 { self.log }
}

/// Prints `d(S³_{p/q}(K), i)` for `0 ≤ i < |p|`, one line each.
pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    let (p, q) = (*args.slope.numer(), *args.slope.denom());
    if p == 0 {
        return Err("0-surgery is not supported.".into());
    }

    let gd = load_diagram(&args.knot)?;
    let gd = if args.mirror { gd.rotate90() } else { gd };
    gd.log_diagram();

    let config = G0CalcConfig { accept_if_unique: true, ..args.config.config() };
    match d_with_config(&gd, p, q, args.g, args.g4, args.tau, config) {
        Ok(ds) => Ok(ds.iter().enumerate().map(|(i, d)| format!("{i}: {d}")).collect::<Vec<_>>().join("\n")),
        Err(e) => {
            error!("Unable to compute d: G0 is not decided.");
            Err(format!("{e:?}").into())
        }
    }
}
