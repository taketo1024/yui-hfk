use clap::ValueEnum;
use yui_hfk::{GridComplex, GCTilde, GCHat, GCMinus};
use yui_homology::DisplayTable;
use crate::app::{AppArgs, utils::load_diagram};

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Knot name (e.g. "3_1") or path to JSON file
    pub knot: String,

    /// Grid homology type to compute
    #[arg(short = 't', long = "type", value_enum, default_value = "hat")]
    pub gh_type: GHType,

    /// Use the mirror of the knot
    #[arg(short, long)]
    pub mirror: bool,

    /// Lower bound for the Alexander grading of generators
    #[arg(short = 'A', long)]
    pub a_bound: Option<isize>,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 { self.log }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum GHType {
    Tilde,
    Hat,
    Minus,
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    let d = load_diagram(&args.knot)?;
    let diagram = if args.mirror { d.rotate90() } else { d };

    let a_bound = args.a_bound;
    let pred = move |x: &yui_hfk::GCGen| {
        a_bound.map_or(true, |b| x.alex_deg >= b)
    };

    let mut buff = String::new();

    match args.gh_type {
        GHType::Tilde => {
            buff.push_str("GH~ (Tilde):\n");
            let h = GridComplex::<GCTilde>::new_filtered(diagram, pred).into_homology();
            buff.push_str(&h.display_table());
        }
        GHType::Hat => {
            buff.push_str("GH^ (Hat):\n");
            let h = GridComplex::<GCHat>::new_filtered(diagram, pred).into_homology();
            buff.push_str(&h.display_table());
        }
        GHType::Minus => {
            buff.push_str("GH- (Minus):\n");
            let h = GridComplex::<GCMinus>::new_filtered(diagram, pred).into_homology();
            buff.push_str(&h.display_table());
        }
    }

    Ok(buff)
}

