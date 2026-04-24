use yui_hfk::GridDiagram;

/// Load a GridDiagram from:
/// - arc-presentation JSON array: `[1,2,3,1,2,3]`
/// - inline sequences JSON: `{"O":[...],"X":[...]}`
/// - knot name: `"3_1"` or file path to a JSON file
pub fn load_diagram(input: &str) -> Result<GridDiagram, Box<dyn std::error::Error>> {
    if let Ok(code) = serde_json::from_str::<Vec<u8>>(input) {
        return Ok(GridDiagram::from_arc_presentation(&code));
    }
    if let Ok(pair) = serde_json::from_str::<[Vec<u8>; 2]>(input) {
        return Ok(GridDiagram::from_sequences(&pair[0], &pair[1]));
    }
    GridDiagram::load(input)
}

#[cfg(test)]
mod tests {
    use super::*;
    use yui_hfk::GridDiagram;

    #[test]
    fn test_load_by_name() {
        let d = load_diagram("3_1").unwrap();
        assert_eq!(d, GridDiagram::trefoil());
    }

    #[test]
    fn test_load_by_path() {
        let dir = env!("CARGO_MANIFEST_DIR");
        let path = format!("{dir}/resources/3_1.json");
        let d = load_diagram(&path).unwrap();
        assert_eq!(d, GridDiagram::trefoil());
    }

    #[test]
    fn test_load_by_arc_presentation() {
        let d = load_diagram("[3,5,2,4,1,3,5,2,4,1]").unwrap();
        assert_eq!(d, GridDiagram::trefoil());
    }

    #[test]
    fn test_load_by_sequences_json() {
        let d = load_diagram("[[2,1,0,4,3],[4,3,2,1,0]]").unwrap();
        assert_eq!(d, GridDiagram::trefoil());
    }

    #[test]
    fn test_load_invalid() {
        assert!(load_diagram("not_a_valid_input").is_err());
    }
}
