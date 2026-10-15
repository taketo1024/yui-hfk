use clap::{Parser, Subcommand};
use crate::app::cmd::{d, tau, upsilon};

use super::args::AppArgs;
use super::cmd::{g0, gh, hfk, setup};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    #[command(subcommand)]
    pub command: Cmd,
}

#[derive(Subcommand, Debug)]
#[clap(rename_all = "lower")]
pub enum Cmd {
    G0(g0::Args),
    Gh(gh::Args),
    Hfk(hfk::Args),
    Tau(tau::Args),
    /// The Υ-invariant, from G0
    Upsilon(upsilon::Args),
    /// The d-invariants of the p/q-surgery, from the V_k
    D(d::Args),
    /// Download grid diagrams from KnotInfo into the data directory
    Setup(setup::Args),
}

impl CliArgs {
    fn app_args(&self) -> &dyn AppArgs {
        match &self.command {
            Cmd::G0(args)  => args,
            Cmd::Gh(args)  => args,
            Cmd::Hfk(args) => args,
            Cmd::Tau(args) => args,
            Cmd::Upsilon(args) => args,
            Cmd::D(args) => args,
            Cmd::Setup(args) => args,
        }
    }

    fn log_level(&self) -> log::LevelFilter { 
        self.app_args().log_level()
    }
}

pub struct App {
    pub args: CliArgs,
}

impl App {
    pub fn new() -> Self {
        App { args: CliArgs::parse() }
    }

    pub fn run(&self) -> Result<String, Box<dyn std::error::Error>> {
        self.init_logger();

        let start = std::time::Instant::now();

        let res = match &self.args.command {
            Cmd::G0(args)  => g0::dispatch(args),
            Cmd::Gh(args)  => gh::dispatch(args),
            Cmd::Hfk(args) => hfk::dispatch(args),
            Cmd::Tau(args) => tau::dispatch(args),
            Cmd::Upsilon(args) => upsilon::dispatch(args),
            Cmd::D(args) => d::dispatch(args),
            Cmd::Setup(args) => setup::dispatch(args),
        };

        log::info!("run time: {}", format_duration(start.elapsed()));

        res
    }

    fn init_logger(&self) {
        let level = self.args.log_level();
        env_logger::Builder::new().filter_level(level).parse_default_env().init();
    }
}

fn format_duration(d: std::time::Duration) -> String {
    let total = d.as_secs();
    let (h, m, s) = (total / 3600, (total % 3600) / 60, total % 60);
    match (h, m) {
        (0, 0) => format!("{:.3}s", d.as_secs_f64()),
        (0, _) => format!("{}m {}s", m, s),
        _      => format!("{}h {}m {}s", h, m, s),
    }
}
