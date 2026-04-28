use clap::{Parser, Subcommand};
use super::args::AppArgs;
use super::cmd::{g0, gh, hfk};

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
}

impl CliArgs {
    fn app_args(&self) -> &dyn AppArgs {
        match &self.command {
            Cmd::G0(args)  => args,
            Cmd::Gh(args)  => args,
            Cmd::Hfk(args) => args,
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
        match &self.args.command {
            Cmd::G0(args)  => g0::dispatch(args),
            Cmd::Gh(args)  => gh::dispatch(args),
            Cmd::Hfk(args) => hfk::dispatch(args),
        }
    }

    fn init_logger(&self) {
        let level = self.args.log_level();
        yui_core::util::log::init_simple_logger(level).unwrap();
    }
}
