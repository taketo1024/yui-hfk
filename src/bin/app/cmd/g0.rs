use yui_hfk::{G0Calc, G0CalcConfig};
use crate::app::{AppArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// Use the mirror of the knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Skip HFK-based region filtering
    #[arg(long)]
    pub no_filter_regions: bool,

    /// Skip generator filtering
    #[arg(long)]
    pub no_filter_gens: bool,

    /// Expected genus; aborts if it disagrees with the computed value
    #[arg(long)]
    pub g: Option<usize>,

    /// Filter candidates using the 4-genus bound
    #[arg(long)]
    pub g4: Option<usize>,

    /// Filter candidates using a known value of the tau invariant
    #[arg(long)]
    pub tau: Option<isize>,

    /// Try the U-specialization NG fast-path on large systems
    #[arg(long)]
    pub try_fast_ng: bool,

    /// Resign on linear systems whose dimension exceeds this size
    #[arg(long)]
    pub max_size: Option<usize>,

    /// Cap on the sparse pre-PLUQ pivot count in the linear solver
    #[arg(long)]
    pub max_piv: Option<usize>,

    /// Schur-complement chunk size in the incremental solver
    #[arg(long)]
    pub chunk: Option<usize>,

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

    let defaults = G0CalcConfig::default();
    let config = G0CalcConfig {
        filter_regions: !args.no_filter_regions,
        filter_gens:    !args.no_filter_gens,
        try_fast_ng:    args.try_fast_ng,
        max_size:       args.max_size.unwrap_or(defaults.max_size),
        max_piv:        args.max_piv.unwrap_or(defaults.max_piv),
        chunk:          args.chunk.unwrap_or(defaults.chunk),
    };

    let mut calc = G0Calc::new(gd, args.g4, args.tau, config);

    if let Some(g) = args.g {
        if g != calc.genus() {
            return Err(format!(
                "expected genus g = {g}, but computed g = {}",
                calc.genus()
            ).into());
        }
    }

    calc.run();

    let format_lines = |rs: &[yui_hfk::g_zero::Region]| -> String {
        rs.iter().map(|r| format!("  {r}")).collect::<Vec<_>>().join("\n")
    };

    let result = if calc.resigned().is_empty() {
        calc.accepted().iter().map(|r| r.to_string()).collect::<Vec<_>>().join("\n")
    } else {
        format!(
            "accepted:\n{}\n\nresigned:\n{}",
            format_lines(calc.accepted()),
            format_lines(calc.resigned()),
        )
    };

    Ok(result)
}
