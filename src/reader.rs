use std::error::Error;
use std::io::BufRead;

use noodles::fasta::{
    record::{Definition, Sequence},
    Record,
};
use regex::Regex;

use crate::msa::MSA;

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn read_clustal(&mut self) -> Result<MSA, Box<dyn Error>> {
        read_clustal(&mut self.inner)
    }
}

fn read_clustal<R>(reader: &mut R) -> Result<MSA, Box<dyn Error>>
where
    R: BufRead,
{
    let mut msa = MSA::default();
    let reader = Reader::new(reader);

    // Parsing program name
    let mut buf = String::new();
    reader.inner.read_line(&mut buf)?;

    let known_header = ["CLUSTAL", "PROBCONS", "MUSCLE", "MSAPROBS", "Kalign"];

    if let Some(header) = known_header.iter().find(|&&h| buf.starts_with(h)) {
        msa.add_annotation("program".to_string(), header.to_string());
    } else {
        return Err(format!(
            "The header is not a recognised CLUSTAL header: {}",
            known_header.join(",")
        )
        .into());
    }

    // Parsing program version
    let version_re = Regex::new(r"(\d+(?:\.\d+)+)").unwrap();
    if let Some(version) = version_re.captures(&buf) {
        msa.add_annotation("version".to_string(), version[1].to_string());
    }

    // Handling rest of file
    let mut consensus_idx: u32 = 1;
    let mut start: usize = 0;
    let mut end: usize = 0;
    buf.clear();
    while reader.inner.read_line(&mut buf)? != 0 {
        if !buf.starts_with(' ') && buf != "\n" {
            let fields: Vec<&str> = buf.split_whitespace().filter(|x| !x.is_empty()).collect();
            start = fields[0].len() + buf[fields[0].len()..].find(fields[1]).unwrap_or(0);
            end = start + fields[1].len();
            let record = Record::new(
                Definition::new(fields[0].to_string(), None),
                Sequence::from(fields[1].as_bytes().to_vec()),
            );
            msa.add_record(record.clone());
        }

        if buf.starts_with(' ') {
            msa.add_column_annotation(consensus_idx.to_string(), buf[start..end].to_string());
        }
        buf.clear();
        consensus_idx += 1;
    }

    Ok(msa)
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::BufReader};

    use super::*;

    #[test]
    fn test_clustal() {
        let mut data = Reader::new(BufReader::new(File::open("tests/clustalw.aln").unwrap()));
        let msa = data.read_clustal().unwrap();
        assert_eq!(msa.get_annotation("program").unwrap(), "CLUSTAL");
        assert_eq!(msa.get_annotation("version").unwrap(), "1.81");
    }
}
