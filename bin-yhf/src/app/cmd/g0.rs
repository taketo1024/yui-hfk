use yui_hfk::{G0Calc, G0CalcConfig};
use yui_hfk::g_zero::{Candidates, Mode, Region};
use crate::app::{AppArgs, utils::load_diagram};

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

    /// Filter candidates using a known value of the tau invariant
    #[arg(short, long)]
    pub tau: Option<isize>,

    /// Test only this region, e.g. "[(0, 1), (1, 0)]"
    #[arg(short='R', long, conflicts_with_all = ["g", "g4", "tau"])]
    pub region: Option<Region>,

    /// A region expected in G0, e.g. "R(0, 1)" (repeatable). Expected candidates are
    /// tested first, searching for a solution first; the results do not depend on it.
    #[arg(long = "expect", value_name = "REGION")]
    pub expect: Vec<Region>,

    /// Skip HFK-based region filtering
    #[arg(long)]
    pub no_filter_regions: bool,

    /// Skip generator filtering
    #[arg(long)]
    pub no_filter_gens: bool,

    /// Don't auto-accept when only one candidate region remains
    #[arg(long)]
    pub no_accept_if_unique: bool,

    /// Try the U-specialization NG fast-path on large systems
    #[arg(long)]
    pub try_fast_ng: bool,

    /// Solver for the linear systems: lanczos or pluq
    #[arg(long, default_value = "lanczos")]
    pub mode: Mode,

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
    validate(args);

    let gd = load_diagram(&args.knot)?;
    let gd = if args.mirror { gd.rotate90() } else { gd };

    let defaults = G0CalcConfig::default();
    let config = G0CalcConfig {
        filter_regions: !args.no_filter_regions,
        filter_gens:    !args.no_filter_gens,
        accept_if_unique: !args.no_accept_if_unique,
        try_fast_ng:    args.try_fast_ng,
        mode:           args.mode,
        max_size:       args.max_size.unwrap_or(defaults.max_size),
        max_piv:        args.max_piv.unwrap_or(defaults.max_piv),
        chunk:          args.chunk.unwrap_or(defaults.chunk),
    };

    let mut calc = match &args.region {
        Some(r) => {
            let mut candidates = Candidates::new();
            candidates.add(r.clone());
            G0Calc::init_with_candidates(&gd, candidates, config)
        }
        None => G0Calc::init(&gd, args.g, args.g4, args.tau, config),
    };

    calc.set_expected(args.expect.clone());
    calc.run();

    if let Some(r) = &args.region {
        return Ok(format!("{r}: {}", region_status(&calc, r)));
    }

    let format_lines = |rs: &[Region]| -> String {
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

// note: a rejected region is dropped without a record, so "neither accepted nor resigned" means NG.
fn region_status(calc: &G0Calc, r: &Region) -> &'static str {
    if calc.accepted().contains(r) {
        "OK"
    } else if calc.resigned().contains(r) {
        "RESIGN"
    } else {
        "NG"
    }
}

fn validate(args: &Args) { 
    if let Some(g) = args.g && let Some(g4) = args.g4 { 
        assert!(g4 <= g, "Invalid inequality: g4 = {g4} > g = {g}");
    }

    if let Some(g4) = args.g4 && let Some(tau) = args.tau { 
        let t_abs = tau.unsigned_abs();
        assert!(t_abs <= g4, "Invalid inequality: |tau| = {t_abs} > g4 = {g4}");
    }
}