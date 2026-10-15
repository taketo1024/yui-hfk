use yui_hfk::G0CalcConfig;
use yui_hfk::g_zero::Mode;

pub trait AppArgs {
    fn log(&self) -> u8;

    fn log_level(&self) -> log::LevelFilter {
        use log::LevelFilter::*;
        match self.log() {
            1 => Info,
            2 => Debug,
            3 => Trace,
            _ => Off,
        }
    }
}

/// The options of `G0CalcConfig`, shared by the commands.
#[derive(Clone, Debug, clap::Args)]
pub struct ConfigArgs {
    /// Try the U-specialization NG fast-path on large systems
    #[arg(long)]
    pub try_fast_ng: bool,

    /// Solver for the linear systems: lanczos or pluq
    #[arg(long, default_value = "lanczos")]
    pub mode: Mode,

    /// Resign on linear systems whose dimension exceeds this size
    #[arg(long)]
    pub max_size: Option<usize>,
}

impl ConfigArgs {
    // note: accept_if_unique stays on, as it always has been in the CLI.
    pub fn config(&self) -> G0CalcConfig {
        let defaults = G0CalcConfig::default();
        G0CalcConfig {
            accept_if_unique: true,
            try_fast_ng: self.try_fast_ng,
            mode:        self.mode,
            max_size:    self.max_size.unwrap_or(defaults.max_size),
            ..defaults
        }
    }
}
