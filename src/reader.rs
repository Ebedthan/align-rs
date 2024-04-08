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

    pub fn read_clustal(&mut self) -> Result<MSA, String> {
        read_clustal(&mut self.inner)
    }
}

fn read_clustal<R>(reader: &mut R) -> Result<MSA, String>
where
    R: BufRead,
{
    let mut msa = MSA::default();
    let reader = Reader::new(reader);

    // Parsing program name
    let mut buf = String::new();
    reader
        .inner
        .read_line(&mut buf)
        .expect("reading first line");

    let known_header = ["CLUSTAL", "PROBCONS", "MUSCLE", "MSAPROBS", "Kalign"];

    if let Some(header) = known_header.into_iter().next() {
        if buf.starts_with(header) {
            msa.add_annotation("program".to_string(), header.to_string());
        } else {
            return Err(format!(
                "{} is not a known CLUSTAL header: {}",
                header,
                known_header.join(",")
            ));
        }
    }

    // Parsing program version
    let version_re = Regex::new(r"(\d+(?:\.\d+)+)").unwrap();
    if let Some(version) = version_re.captures(&buf) {
        msa.add_annotation("version".to_string(), version[1].to_string());
    }

    // handling the two empty lines
    buf.clear();
    reader.inner.read_line(&mut buf).expect("read line");
    while buf.trim_start_matches(' ').trim_end_matches(' ') == "\n" {
        reader.inner.read_line(&mut buf).expect("read line");
    }

    // Handling rest of file
    let mut consensus_idx = 1;
    loop {
        buf.clear();
        reader.inner.read_line(&mut buf).expect("read line");
        let mut start = 0;
        let mut end = 0;
        if buf.chars().nth(0).unwrap() != ' '
            && buf.trim_start_matches('\n').trim_end_matches('\n') != ""
        {
            println!("1");
            let fields = buf.trim_end_matches(' ').split(' ').collect::<Vec<&str>>();
            println!("{:?}", fields);
            start = fields[0].len() + buf[fields[0].len()..].find(fields[1]).unwrap_or(0);
            end = start + fields[1].len();
            let record = Record::new(
                Definition::new(fields[0].to_string(), None),
                Sequence::from(fields[1].as_bytes().to_vec()),
            );
            msa.add_record(record.clone()).expect("sequence should not already exists and should have the same length as existing sequences");
        } else if buf.chars().nth(0).unwrap() != ' ' {
            println!("2");
            msa.add_column_annotation(consensus_idx.to_string(), buf[start..end].to_string());
        } else {
            println!("3");
            break;
        }
        consensus_idx += 1;
        println!("{}", consensus_idx);
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
        println!(
            "program: {}\nversion: {}",
            msa.get_annotation("program").unwrap(),
            msa.get_annotation("version").unwrap()
        );
        println!("{:?}", msa.get_ids());
    }
}
