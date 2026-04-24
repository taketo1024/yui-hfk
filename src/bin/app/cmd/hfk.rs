use yui_hfk::HFK;
use crate::app::{AppArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// Use the mirror of each knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 { self.log }
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    let gd= load_diagram(&args.knot)?;
    let gd = if args.mirror { gd.rotate90() } else { gd };
    let h = HFK::from(&gd);
    Ok(h.display_table("M", "A"))
}
