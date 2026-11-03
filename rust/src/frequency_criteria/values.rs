//! Frequency-owned typed value boundary. Host ownership and mature dependencies remain explicit.
pub type Result<T> = std::result::Result<T, String>;

#[derive(Clone, Debug, PartialEq)]
pub enum Cell {
    None,
    Na,
    Bool(bool),
    Int(i64),
    Float(f64),
    Text(String),
}

impl Cell {
    pub(crate) fn missing(&self) -> bool {
        matches!(self, Self::None | Self::Na)
            || matches!(self, Self::Float(value) if value.is_nan())
    }

    pub(crate) fn number(&self) -> Option<f64> {
        match self {
            Self::Bool(value) => Some(f64::from(u8::from(*value))),
            Self::Int(value) => Some(*value as f64),
            Self::Float(value) => Some(*value),
            _ => None,
        }
    }

    pub(crate) fn truth(&self) -> Result<Option<bool>> {
        Ok(match self {
            Self::Na => None,
            Self::None => Some(false),
            Self::Bool(value) => Some(*value),
            Self::Int(value) => Some(*value != 0),
            Self::Float(value) => Some(*value != 0.0),
            Self::Text(value) => Some(!value.is_empty()),
        })
    }
}

/// A typed whole column, ndarray, list, or scalar. Axis keys are lossless
/// host-encoded ordinary labels, not numeric positions masquerading as labels.
#[derive(Clone, Debug, PartialEq)]
pub struct Value {
    pub kind: String,
    pub dtype: String,
    pub index: Vec<String>,
    pub cells: Vec<Cell>,
}

impl Value {
    pub fn scalar(cell: Cell) -> Self {
        Self {
            kind: "scalar".into(),
            dtype: "scalar".into(),
            index: Vec::new(),
            cells: vec![cell],
        }
    }

    pub(crate) fn column(&self, cells: Vec<Cell>, dtype: &str) -> Result<Self> {
        if self.kind != "scalar" && self.cells.len() != cells.len() {
            return Err("column result length differs".into());
        }
        Ok(Self {
            kind: self.kind.clone(),
            dtype: dtype.into(),
            index: self.index.clone(),
            cells,
        })
    }

    pub(crate) fn text(&self) -> Result<&str> {
        match self.cells.as_slice() {
            [Cell::Text(value)] if self.kind == "scalar" => Ok(value),
            _ => Err("expected scalar text".into()),
        }
    }
}

