use knuckles_parse::records::Record;
use std::collections::HashMap;
use std::sync::OnceLock;
use wasm_bindgen::prelude::*;

fn terms_map() -> &'static HashMap<&'static str, &'static str> {
    static TERMS: OnceLock<HashMap<&'static str, &'static str>> = OnceLock::new();
    TERMS.get_or_init(|| {
        let mut m = HashMap::new();
        m.insert("ALA", "A");
        m.insert("VAL", "V");
        m.insert("PHE", "F");
        m.insert("PRO", "P");
        m.insert("MET", "M");
        m.insert("ILE", "I");
        m.insert("LEU", "L");
        m.insert("ASP", "D");
        m.insert("GLU", "E");
        m.insert("LYS", "K");
        m.insert("ARG", "R");
        m.insert("SER", "S");
        m.insert("THR", "T");
        m.insert("TYR", "Y");
        m.insert("HIS", "H");
        m.insert("CYS", "C");
        m.insert("ASN", "N");
        m.insert("GLN", "Q");
        m.insert("TRP", "W");
        m.insert("GLY", "G");
        m
    })
}

fn aa3to1(aa: &str) -> &'static str {
    terms_map().get(aa).copied().unwrap_or("X")
}

fn handle_record(
    writer: &mut impl std::io::Write,
    record: Record,
    prev_chain_id: &mut char,
) -> Result<(), std::io::Error> {
    match record {
        Record::Atom(record) => {
            let name = record.name;
            let alt_loc = record.alt_loc;
            let chain_id = record.chain_id.unwrap();
            if alt_loc.is_some() && alt_loc != Some('A') || !name.is_empty() && name != "CA" {
                return Ok(());
            }
            if *prev_chain_id != chain_id {
                *prev_chain_id = chain_id;
                writeln!(writer, ">pdb:{}", chain_id)?;
            }
            let res_name: &str = &record.res_name;
            let char = aa3to1(res_name);
            if char == "X" {
                eprintln!("Unknown amino acid: {}", res_name);
            }
            write!(writer, "{}", char)?;
            Ok(())
        }
        Record::Hetatm(record) => {
            let name = record.name;
            let alt_loc = record.alt_loc;
            let res_name = record.res_name;
            let chain_id = record.chain_id.unwrap();
            if alt_loc.is_some() && alt_loc != Some('A')
                || !name.is_empty() && name != "CA"
                || !res_name.is_empty() && res_name != "MSE"
            {
                return Ok(());
            }
            if *prev_chain_id != chain_id {
                *prev_chain_id = chain_id;
                writeln!(writer, ">pdb:{}", chain_id)?;
            }
            let char = aa3to1(&res_name);
            write!(writer, "{}", char)?;
            Ok(())
        }
        _ => Ok(()),
    }
}

fn pdb_to_fasta(writer: &mut impl std::io::Write, contents: String) -> Result<(), std::io::Error> {
    let mut prev_chain_id = ' ';
    for line in contents.lines() {
        if line.len() < 6 {
            continue;
        }
        let record = Record::try_from(line);
        if record.is_err() {
            continue;
        }
        let record = record.unwrap();
        match record {
            Record::Endmdl() => break,
            Record::Atom(_) | Record::Hetatm(_) => {
                handle_record(writer, record, &mut prev_chain_id)?;
            }
            _ => {}
        }
    }

    writer.flush()
}

#[wasm_bindgen]
pub fn pdb_to_fasta_js(contents: String) -> String {
    let mut out = Vec::new();
    let _ = pdb_to_fasta(&mut out, contents);
    String::from_utf8(out).unwrap()
}

// fn main() {
//     let mut out = std::io::stdout();
//     let contents = std::fs::read_to_string("4pth.pdb").unwrap();
//     let _ = pdb_to_fasta(&mut out, contents);
//     // println!("Hello, world!");
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aa3to1() {
        let input = "ALA";
        assert_eq!(aa3to1(input), "A");
        let input = "OTO";
        assert_eq!(aa3to1(input), "X");
    }
}
