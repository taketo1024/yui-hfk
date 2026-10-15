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

    /// Solver for the linear systems: lanczos, lanczos(rounds), pluq or pluq(chunk)
    #[arg(long, default_value = "lanczos", value_parser = parse_mode)]
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

/// `lanczos` and `pluq` with the default values, or `lanczos(n)` and `pluq(n)` with `n > 0`.
fn parse_mode(s: &str) -> Result<Mode, String> {
    let err = || format!("expected `lanczos`, `lanczos(n)`, `pluq` or `pluq(n)` with n > 0, got `{s}`");
    let (name, n) = match s.strip_suffix(')').and_then(|t| t.split_once('(')) {
        Some((name, n)) => {
            let n = n.replace('_', "").parse::<usize>().ok().filter(|&n| n > 0).ok_or_else(err)?;
            (name, Some(n))
        }
        None => (s, None),
    };
    match (name, n) {
        ("lanczos", None) => Ok(Mode::default()),
        ("lanczos", Some(rounds)) => Ok(Mode::lanczos(rounds)),
        ("pluq", None) => Ok(Mode::pluq(100_000)),
        ("pluq", Some(chunk)) => Ok(Mode::pluq(chunk)),
        _ => Err(err()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_mode() {
        assert_eq!(parse_mode("lanczos"), Ok(Mode::lanczos(1)));
        assert_eq!(parse_mode("lanczos(3)"), Ok(Mode::lanczos(3)));
        assert_eq!(parse_mode("pluq"), Ok(Mode::pluq(100_000)));
        assert_eq!(parse_mode("pluq(50_000)"), Ok(Mode::pluq(50_000)));
        for s in ["foo", "lanczos()", "lanczos(0)", "pluq(x)", "pluq(1", "lanczos 3"] {
            assert!(parse_mode(s).is_err(), "{s}");
        }
    }
}
