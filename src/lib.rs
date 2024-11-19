use knuckles_parse::records::Record;
use wasm_bindgen::prelude::*;

static TERMS: phf::Map<&'static str, &'static str> = phf::phf_map! {
    "ALA" => "A",
    "VAL" => "V",
    "PHE" => "F",
    "PRO" => "P",
    "MET" => "M",
    "ILE" => "I",
    "LEU" => "L",
    "ASP" => "D",
    "GLU" => "E",
    "LYS" => "K",
    "ARG" => "R",
    "SER" => "S",
    "THR" => "T",
    "TYR" => "Y",
    "HIS" => "H",
    "CYS" => "C",
    "ASN" => "N",
    "GLN" => "Q",
    "TRP" => "W",
    "GLY" => "G",
};

fn aa3to1(aa: &str) -> &'static str {
    TERMS.get(aa).unwrap_or(&"X")
    // terms_map().get(aa).copied().unwrap_or("X")
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
