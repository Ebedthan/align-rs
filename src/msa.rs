use core::fmt;
use std::collections::HashMap;

use crate::record::Record;

/// Structure containing multiple sequence alignments
///
#[derive(Default, Debug, Clone, PartialEq)]
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

    pub fn col_len(&self) -> usize {
        self.records.first().map(|x| x.len()).unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn clear(&mut self) {
        self.records.clear();
        self.annotations.clear();
        self.column_annotations.clear();
    }

    pub fn get_annotation(&self, name: &str) -> Option<&String> {
        self.annotations.get(name)
    }

    pub fn get_column_annotation(&self, name: &str) -> Option<&String> {
        self.column_annotations.get(name)
    }

    pub fn add_column_annotation(&mut self, name: &str, value: &str) {
        self.column_annotations
            .entry(name.to_string())
            .or_default()
            .push_str(value);
    }

    pub fn add_annotation(&mut self, name: String, value: String) -> Option<String> {
        self.annotations.insert(name, value)
    }

    pub fn contains(&self, haystack: &str) -> bool {
        self.records.iter().map(|x| x.id()).any(|x| x == haystack)
    }

    pub fn push_record(&mut self, id: &str, seq: &str) {
        if self.contains(id) {
            for x in &mut self.records {
                if x.id() == id {
                    x.push_seq(seq);
                    break;
                }
            }
        } else {
            self.records.push(Record::new(id, seq));
        }
    }
}

impl fmt::Display for MSA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_empty() {
            return write!(f, "No sequence in alignment");
        }

        let num_rows = self.len();
        let num_cols = self.col_len();
        let mut string = String::new();
        let mut to_continue: &str = "";
        for (idx, record) in self.records.iter().enumerate() {
            if idx > 9 {
                to_continue = "...";
                break;
            }
            let truncated_sequence = &record.sequence()[..std::cmp::min(30, record.len())];
            string.push_str(&format!(
                "{}\t{}\n",
                record.id(),
                if record.len() > 30 {
                    format!("{}...", truncated_sequence)
                } else {
                    truncated_sequence.to_string()
                }
            ));
        }
        let header = if num_rows == 1 {
            format!("Alignment with {} row", num_rows)
        } else {
            format!("Alignment with {} rows", num_rows)
        };

        let column_desc = if num_cols == 1 {
            format!(" and {} column", num_cols)
        } else {
            format!(" and {} columns", num_cols)
        };
        if to_continue.is_empty() {
            write!(f, "{}{}\n{}", header, column_desc, string)
        } else {
            write!(f, "{}{}\n{}\n{}", header, column_desc, string, to_continue)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Record;

    #[test]
    fn msa_new() {
        let records = vec![Record::new("id1", "ACGT"), Record::new("id2", "TGCA")];

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
        let records = vec![Record::new("id1", "ACGT"), Record::new("id2", "TGCA")];

        let msa = MSA::new(records.clone(), HashMap::new(), HashMap::new());
        assert_eq!(msa.len(), 2);
    }

    #[test]
    fn msa_alignment_len() {
        let records = vec![Record::new("id1", "ACGT"), Record::new("id2", "TGCA")];

        let msa = MSA::new(records.clone(), HashMap::new(), HashMap::new());
        assert_eq!(msa.col_len(), 4);
    }

    #[test]
    fn msa_push_valid_sequence() {
        let mut msa = MSA::default();
        msa.push_record("id1", "ACGT");
        assert_eq!(msa.len(), 1);
        assert_eq!(msa.records[0], Record::new("id1", "ACGT"));
    }

    #[test]
    fn msa_print_no_seqs() {
        let msa = MSA::default();
        assert_eq!(msa.to_string(), "No sequence in alignment");
    }
    #[test]
    fn msa_print_one_seqs() {
        let mut msa = MSA::default();
        msa.push_record("id1", "ACG");
        assert_eq!(
            msa.to_string(),
            "Alignment with 1 row and 3 columns\nid1\tACG\n"
        );
    }
}
