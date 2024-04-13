use std::collections::HashMap;

/// Simple Sequence Record Structure for multiple sequence alignment
///
#[derive(Debug, Clone, Default)]
pub struct Record {
    /// Sequence ID
    id: String,

    /// Sequence string
    sequence: String,

    /// Letter annotation
    annotation: HashMap<String, String>,
}

impl Record {
    /// Return sequence length
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    pub fn new(id: &str, sequence: &str) -> Record {
        Record {
            id: id.to_string(),
            sequence: sequence.to_string(),
            annotation: HashMap::new(),
        }
    }

    pub fn push_seq(&mut self, string: &str) {
        self.sequence.push_str(string);
    }

    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.sequence.is_empty()
    }

    pub fn push_record(&mut self, id: &str, seq: &str) {
        if self.is_empty() {
            self.id = id.to_string();
            self.sequence = seq.to_string();
        } else {
            self.sequence.push_str(seq);
        }
    }
}
