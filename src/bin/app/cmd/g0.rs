use yui_hfk::G0Calc;
use crate::app::{AppArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// Use the mirror of the knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Filter candidates using the 4-genus bound
    #[arg(long)]
    pub g4: Option<usize>,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 { self.log }
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    let gd = load_diagram(&args.knot)?;
    let gd = if args.mirror { gd.rotate90() } else { gd };

    let mut calc = G0Calc::new(gd, true, true);
    if let Some(g4) = args.g4 {
        calc.filter_by_g4(g4);
    }
    calc.run();

    let result = calc.accepted().iter()
        .map(|r| r.to_string())
        .collect::<Vec<_>>()
        .join("\n");

    Ok(result)
}
