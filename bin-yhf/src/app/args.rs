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
