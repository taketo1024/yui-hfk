use yhf::App;

fn main() {
    let app = App::new();
    match app.run() {
        Ok(output) => println!("{output}"),
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    }
}
