use core::fmt;
use std::collections::HashMap;

use noodles::fasta::Record;

#[derive(Default, Debug, Clone)]
pub struct MSA {
    /// A list of Record objects, whose sequences are all the same length
    records: Vec<Record>,

    /// Information about the whole alignment
    annotations: HashMap<String, String>,

    /// Per column annotation
    column_annotations: HashMap<String, String>,
}

impl MSA {
    /// Creates a new multiple sequence alignment structure
    pub fn new(
        records: Vec<Record>,
        annotations: HashMap<String, String>,
        column_annotations: HashMap<String, String>,
    ) -> Self {
        MSA {
            records,
            annotations,
            column_annotations,
        }
    }

    /// Returns the number of record in alignment
    /// # Example
    /// ```
    /// use align_rs::msa::MSA;
    ///
    /// let msa = MSA::default();
    /// assert_eq!(msa.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn alignment_len(&self) -> usize {
        if self.records.is_empty() {
            0
        } else {
            self.records[0].sequence().len()
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn get_ids(&self) -> Vec<String> {
        self.records
            .iter()
            .map(|x| String::from_utf8(x.name().to_vec()).unwrap())
            .collect()
    }

    pub fn add_record(&mut self, record: Record) {
        self.records.push(record);
    }

    pub fn clear(&mut self) {
        self.records.clear();
        self.annotations.clear();
        self.column_annotations.clear();
    }

    pub fn records(&self) -> &Vec<noodles::fasta::Record> {
        &self.records
    }

    pub fn get_annotation(&self, name: &str) -> Option<&String> {
        self.annotations.get(name)
    }

    pub fn get_column_annotations(&self) -> &HashMap<String, String> {
        &self.column_annotations
    }

    pub fn add_column_annotation(&mut self, name: String, value: String) -> Option<String> {
        self.column_annotations.insert(name, value)
    }

    pub fn add_annotation(&mut self, name: String, value: String) -> Option<String> {
        self.annotations.insert(name, value)
    }
}

impl fmt::Display for MSA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut string: Vec<String> = vec![];
        for record in &self.records {
            let mstr: String = if record.sequence().len() > 60 {
                String::from_utf8(record.name().to_vec()).unwrap()
                    + &String::from(" ")
                    + &String::from_utf8(record.sequence().as_ref()[..60].to_vec()).unwrap()
                    + &String::from("...")
            } else {
                String::from_utf8(record.name().to_vec()).unwrap()
                    + &String::from(" ")
                    + &String::from_utf8(record.sequence().as_ref().to_vec()).unwrap()
            };
            string.push(mstr);
        }
        let mstr = string
            .iter()
            .map(|x| x.to_string() + "\n")
            .collect::<String>();
        let mstr = mstr.trim_end_matches('\n');
        if self.records.len() == 1 && self.records[0].sequence().len() == 1 {
            write!(
                f,
                "Alignment with {} row and {} column\n{}",
                self.records.len(),
                self.records[0].sequence().len(),
                mstr
            )
        } else if self.records.len() == 1 && self.records[0].sequence().len() > 1 {
            write!(
                f,
                "Alignment with {} row and {} columns\n{}",
                self.records.len(),
                self.records[0].sequence().len(),
                mstr
            )
        } else if self.records.is_empty() {
            write!(f, "No sequence in alignment")
        } else {
            write!(
                f,
                "Alignment with {} rows and {} columns\n{}",
                self.records.len(),
                self.records[0].sequence().len(),
                mstr
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::fasta;

    #[test]
    fn msa_new() {
        let records = vec![
            fasta::Record::new(
                fasta::record::Definition::new("id1", None),
                fasta::record::Sequence::from(b"ACGT".to_vec()),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("id2", None),
                fasta::record::Sequence::from(b"TGCA".to_vec()),
            ),
        ];

        let mut annotations = HashMap::new();
        annotations.insert(String::from("author"), String::from("John Doe"));

        let mut column_annotations = HashMap::new();
        column_annotations.insert(String::from("1"), String::from("column 1"));

        let msa = MSA::new(
            records.clone(),
            annotations.clone(),
            column_annotations.clone(),
        );

        assert_eq!(msa.records, records);
        assert_eq!(msa.annotations, annotations);
        assert_eq!(msa.column_annotations, column_annotations);
    }

    #[test]
    fn msa_len() {
        let records = vec![
            fasta::Record::new(
                fasta::record::Definition::new("id1", None),
                fasta::record::Sequence::from(b"ACGT".to_vec()),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("id2", None),
                fasta::record::Sequence::from(b"TGCA".to_vec()),
            ),
        ];

        let msa = MSA::new(records.clone(), HashMap::new(), HashMap::new());
        assert_eq!(msa.len(), 2);
    }

    #[test]
    fn msa_alignment_len() {
        let records = vec![
            fasta::Record::new(
                fasta::record::Definition::new("id1", None),
                fasta::record::Sequence::from(b"ACGT".to_vec()),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("id2", None),
                fasta::record::Sequence::from(b"TGCA".to_vec()),
            ),
        ];

        let msa = MSA::new(records.clone(), HashMap::new(), HashMap::new());
        assert_eq!(msa.alignment_len(), 4);
    }

    #[test]
    fn msa_push_valid_sequence() {
        let mut msa = MSA::default();
        let record = fasta::Record::new(
            fasta::record::Definition::new("id1", None),
            fasta::record::Sequence::from(b"ACGT".to_vec()),
        );
        msa.add_record(record.clone());
        assert_eq!(msa.len(), 1);
        assert_eq!(msa.records[0], record);
    }

    #[test]
    fn msa_print_no_seqs() {
        let msa = MSA::default();
        assert_eq!(msa.to_string(), "No sequence in alignment");
    }
    #[test]
    fn msa_print_one_seqs() {
        let mut msa = MSA::default();
        let record = fasta::Record::new(
            fasta::record::Definition::new("id1", None),
            fasta::record::Sequence::from(b"ACG".to_vec()),
        );
        msa.add_record(record.clone());
        assert_eq!(
            msa.to_string(),
            "Alignment with 1 row and 3 columns\nid1 ACG"
        );
    }
}
