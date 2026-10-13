use std::error::Error;
use std::fs;
use std::io::{Cursor, Read};
use std::path::{Path, PathBuf};

use calamine::{Reader, Xls, open_workbook_from_rs};
use yui_core::util::data_dir::resolve_data_dir;
use zip::ZipArchive;

use crate::app::AppArgs;

const URL: &str = "https://raw.githubusercontent.com/knotinfo-org/assets/main/knotinfo_data_complete.xls.zip";
const XLS_NAME: &str = "knotinfo_data_complete.xls";
const MAX_DOWNLOAD: u64 = 256 << 20;

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Output directory (default: $YUI_DATA_DIR, else the platform user-data dir)
    #[arg(long)]
    pub out: Option<PathBuf>,

    /// Read this local file (.xls or .xls.zip) instead of downloading
    #[arg(long = "in", value_name = "FILE")]
    pub input: Option<PathBuf>,

    /// Overwrite existing files
    #[arg(long)]
    pub force: bool,

    /// Log level: 0=off, 1=info, 2=debug, 3=trace
    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn log(&self) -> u8 {
        self.log
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Outcome {
    Written,
    Kept,
    Skipped,
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn Error>> {
    let out_dir = match &args.out {
        Some(dir) => dir.clone(),
        None => resolve_data_dir()?,
    };
    let xls = match &args.input {
        Some(path) => read_local(path)?,
        None => unzip_xls(&download()?)?,
    };

    let grid_dir = out_dir.join("grid_diagrams");
    fs::create_dir_all(&grid_dir)?;

    eprintln!("==> converting grid notations");
    let outcomes = read_grid_notations(xls)?.iter().map(|(name, grid_notation)|
        save(&grid_dir, name, grid_notation, args.force)
    ).collect::<Result<Vec<_>, _>>()?;

    let count = |o: Outcome| outcomes.iter().filter(|&&x| x == o).count();
    Ok(format!(
        "wrote {} grid diagrams to {} (kept {} existing, skipped {})",
        count(Outcome::Written),
        grid_dir.display(),
        count(Outcome::Kept),
        count(Outcome::Skipped),
    ))
}

fn download() -> Result<Vec<u8>, Box<dyn Error>> {
    eprintln!("==> downloading {URL}");
    let bytes = ureq::get(URL).call()?
        .into_body()
        .into_with_config()
        .limit(MAX_DOWNLOAD)
        .read_to_vec()?;
    Ok(bytes)
}

fn read_local(path: &Path) -> Result<Vec<u8>, Box<dyn Error>> {
    eprintln!("==> reading {}", path.display());
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if path.extension().is_some_and(|ext| ext == "zip") {
        unzip_xls(&bytes)
    } else {
        Ok(bytes)
    }
}

fn unzip_xls(zip: &[u8]) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut archive = ZipArchive::new(Cursor::new(zip))?;
    let mut file = archive.by_name(XLS_NAME)?;
    let mut xls = Vec::with_capacity(file.size() as usize);
    file.read_to_end(&mut xls)?;
    Ok(xls)
}

/// The `(name, grid_notation)` pairs of all rows in the KnotInfo table.
fn read_grid_notations(xls: Vec<u8>) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let mut workbook: Xls<_> = open_workbook_from_rs(Cursor::new(xls))?;
    let sheet = workbook.worksheet_range_at(0).ok_or("the workbook has no sheet")??;

    let mut rows = sheet.rows();
    let header = rows.next().ok_or("the sheet is empty")?;
    let col = |key: &str| {
        header.iter().position(|h| h == key).ok_or_else(|| format!("column `{key}` not found"))
    };
    let i_name = col("name")?;
    let i_grid = col("grid_notation")?;

    // note: the row after the header holds human-readable column titles
    let pairs = rows.skip(1).map(|row| (
        row[i_name].to_string().trim().to_string(),
        row[i_grid].to_string().trim().to_string(),
    )).collect();

    Ok(pairs)
}

fn save(dir: &Path, name: &str, grid_notation: &str, force: bool) -> Result<Outcome, Box<dyn Error>> {
    if name.is_empty() || grid_notation.is_empty() {
        return Ok(Outcome::Skipped);
    }
    if name.contains(['/', '\\']) {
        eprintln!("  skip {name}: invalid name");
        return Ok(Outcome::Skipped);
    }

    let ox = match grid_notation_to_ox(grid_notation) {
        Ok(ox) => ox,
        Err(e) => {
            eprintln!("  skip {name}: {e}");
            return Ok(Outcome::Skipped);
        }
    };

    let path = dir.join(format!("{name}.json"));
    if path.exists() && !force {
        return Ok(Outcome::Kept);
    }
    fs::write(path, serde_json::to_string(&ox)?)?;
    Ok(Outcome::Written)
}

/// Converts KnotInfo's `grid_notation` (1-indexed `[col, row]` markers, two per column)
/// into OX-form `[os, xs]`: `os[c]` (resp. `xs[c]`) is the 0-indexed row of the O (resp. X) in column `c`.
fn grid_notation_to_ox(s: &str) -> Result<[Vec<usize>; 2], String> {
    let points: Vec<[usize; 2]> = serde_json::from_str(s).map_err(|e| e.to_string())?;
    if points.is_empty() || !points.len().is_multiple_of(2) {
        return Err(format!("expected a positive even number of markers, got {}", points.len()));
    }
    let n = points.len() / 2;

    // cols[c] = the two rows with a marker in column c
    let cols = points.chunks(2).enumerate().map(|(c, pair)| {
        let ([c1, r1], [c2, r2]) = (pair[0], pair[1]);
        if c1 != c + 1 || c2 != c + 1 {
            return Err(format!("markers {} and {} must lie in column {}", 2 * c + 1, 2 * c + 2, c + 1));
        }
        if !(1..=n).contains(&r1) || !(1..=n).contains(&r2) || r1 == r2 {
            return Err(format!("invalid rows {r1}, {r2} in column {}", c + 1));
        }
        // note: KnotInfo counts rows in the opposite direction
        Ok([n - r1, n - r2])
    }).collect::<Result<Vec<_>, _>>()?;

    // rows[r] = the two columns with a marker in row r, in increasing order
    let rows = (0..n).map(|r| {
        let cs: Vec<_> = (0..n).filter(|&c| cols[c].contains(&r)).collect();
        <[usize; 2]>::try_from(cs).map_err(|cs| format!("row {} has {} markers, expected 2", n - r, cs.len()))
    }).collect::<Result<Vec<_>, _>>()?;

    let os = assign_os(&cols, &rows);
    let xs = (0..n).map(|c| other(cols[c], os[c])).collect();

    Ok([os, xs])
}

/// The source does not say which marker is the O. For each component, put the O on the left
/// marker of its first row, then walk the component: O → X in the same column → O in the same row.
fn assign_os(cols: &[[usize; 2]], rows: &[[usize; 2]]) -> Vec<usize> {
    let n = cols.len();
    let mut o_of_col = vec![None; n];
    let mut o_of_row = vec![None; n];

    for start in 0..n {
        if o_of_row[start].is_some() {
            continue;
        }
        let (mut c, mut r) = (rows[start][0], start);
        while o_of_row[r].is_none() {
            o_of_col[c] = Some(r);
            o_of_row[r] = Some(c);
            r = other(cols[c], r);
            c = other(rows[r], c);
        }
        debug_assert_eq!(o_of_row[r], Some(c));
    }

    o_of_col.into_iter().map(Option::unwrap).collect()
}

fn other([a, b]: [usize; 2], v: usize) -> usize {
    if v == a { b } else { a }
}

#[cfg(test)]
mod tests {
    use super::*;
    use yui_hfk::GridDiagram;

    #[test]
    fn test_trefoil() {
        let s = "[[1,1],[1,3],[2,2],[2,4],[3,3],[3,5],[4,1],[4,4],[5,2],[5,5]]";
        let [os, xs] = grid_notation_to_ox(s).unwrap();
        assert_eq!(os, vec![2, 1, 0, 4, 3]);
        assert_eq!(xs, vec![4, 3, 2, 1, 0]);

        let to_u8 = |v: &[usize]| v.iter().map(|&i| i as u8).collect::<Vec<_>>();
        assert_eq!(GridDiagram::from_sequences(&to_u8(&os), &to_u8(&xs)), GridDiagram::trefoil());
    }

    #[test]
    fn test_invalid() {
        assert!(grid_notation_to_ox("").is_err());
        assert!(grid_notation_to_ox("[]").is_err());
        assert!(grid_notation_to_ox("[[1,1],[1,2],[2,1]]").is_err());         // odd number of markers
        assert!(grid_notation_to_ox("[[1,1],[2,2],[2,1],[1,2]]").is_err());   // pair not in one column
        assert!(grid_notation_to_ox("[[1,1],[1,1],[2,2],[2,2]]").is_err());   // two markers at one point
        assert!(grid_notation_to_ox("[[1,1],[1,2],[2,1],[2,3]]").is_err());   // row out of range
        assert!(grid_notation_to_ox("[[1,1],[1,2],[2,1],[2,2],[3,1],[3,3]]").is_err()); // row 1 has 3 markers
    }
}
