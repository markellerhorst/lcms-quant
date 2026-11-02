//! .quant project file: plain tab-separated text. Raw data is NOT stored (it is re-imported from the
//! recorded paths and verified by SHA-256); only targets, parameters and per-file analysis state are.
use crate::quant::Analyte;

#[derive(Clone, Debug, PartialEq)]
pub struct CellRec { pub analyte: usize, pub found: bool, pub manual: bool, pub start: f64, pub end: f64, pub include: bool }

#[derive(Clone, Debug, PartialEq)]
pub struct SampleRec {
    pub path: String, pub raw_sha256: String, pub kind: String, pub conc: f64, pub rt_shift: f64, pub show: bool,
    pub cells: Vec<CellRec>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct ProjectFile {
    pub app_version: String,
    pub user: String,
    pub params: Vec<(String, String)>,
    pub analytes: Vec<Analyte>,
    pub samples: Vec<SampleRec>,
}

const MAGIC: &str = "LCMSQUANT-PROJECT";

fn clean(s: &str) -> String { s.replace(['\t', '\n', '\r'], " ") }

impl ProjectFile {
    pub fn to_text(&self) -> String {
        let mut o = format!("{}\t1\n", MAGIC);
        o += &format!("V\t{}\t{}\n", clean(&self.app_version), clean(&self.user));
        for (k, v) in &self.params { o += &format!("P\t{}\t{}\n", clean(k), clean(v)); }
        for a in &self.analytes {
            let rt = a.rt.map(|v| v.to_string()).unwrap_or_else(|| "n/a".into());
            let ions = a.ions.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(";");
            let prec = a.precursor.map(|v| v.to_string()).unwrap_or_else(|| "n/a".into());
            o += &format!("A\t{}\t{}\t{}\t{}\n", clean(&a.name), rt, ions, prec);
        }
        for s in &self.samples {
            o += &format!("S\t{}\t{}\t{}\t{}\t{}\t{}\n", clean(&s.path), s.raw_sha256, clean(&s.kind), s.conc, s.rt_shift, s.show);
            for c in &s.cells {
                o += &format!("C\t{}\t{}\t{}\t{}\t{}\t{}\n", c.analyte, c.found, c.manual, c.start, c.end, c.include);
            }
        }
        o
    }

    pub fn parse(text: &str) -> Result<ProjectFile, String> {
        let mut lines = text.lines();
        let first = lines.next().unwrap_or("");
        if !first.starts_with(MAGIC) { return Err("not an LC-MS/MS Quant project file".into()); }
        let mut p = ProjectFile::default();
        for (n, l) in lines.enumerate() {
            let f: Vec<&str> = l.split('\t').collect();
            let err = |m: &str| format!("line {}: {}", n + 2, m);
            match f[0] {
                "V" if f.len() >= 3 => { p.app_version = f[1].into(); p.user = f[2].into(); }
                "P" if f.len() >= 3 => p.params.push((f[1].into(), f[2].into())),
                "A" if f.len() >= 4 => {
                    let rt = if f[2] == "n/a" { None } else { Some(f[2].parse().map_err(|_| err("bad RT"))?) };
                    let ions: Result<Vec<f64>, _> = f[3].split(';').filter(|s| !s.is_empty()).map(|s| s.parse::<f64>()).collect();
                    let ions = ions.map_err(|_| err("bad ion"))?;
                    if ions.is_empty() { return Err(err("target without ions")); }
                    let precursor = match f.get(4) { Some(s) if *s != "n/a" && !s.is_empty() => Some(s.parse().map_err(|_| err("bad precursor"))?), _ => None };
                    p.analytes.push(Analyte { name: f[1].into(), rt, ions, precursor });
                }
                "S" if f.len() >= 7 => p.samples.push(SampleRec {
                    path: f[1].into(), raw_sha256: f[2].into(), kind: f[3].into(),
                    conc: f[4].parse().map_err(|_| err("bad conc"))?, rt_shift: f[5].parse().map_err(|_| err("bad shift"))?,
                    show: f[6] == "true", cells: vec![],
                }),
                "C" if f.len() >= 7 => {
                    let s = p.samples.last_mut().ok_or_else(|| err("cell before sample"))?;
                    s.cells.push(CellRec {
                        analyte: f[1].parse().map_err(|_| err("bad analyte index"))?, found: f[2] == "true", manual: f[3] == "true",
                        start: f[4].parse().map_err(|_| err("bad start"))?, end: f[5].parse().map_err(|_| err("bad end"))?,
                        include: f[6] == "true",
                    });
                }
                _ => {} // unknown / empty lines are ignored (forward compatible)
            }
        }
        Ok(p)
    }

    pub fn param(&self, key: &str) -> Option<&str> {
        self.params.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }
}
