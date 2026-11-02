use crate::{
    audit::{self, iso_now, Audit},
    fmt::{parse_sci, sci, sci_prec},
    history::History,
    mzml::{self, MzData, Spectrum},
    project::{CellRec, ProjectFile, SampleRec},
    quant::*,
    rawconv,
};
use eframe::egui::{self, Color32, Key, Modifiers, PointerButton, RichText};
use egui_plot::{Legend, Line, Plot, PlotPoint, PlotPoints, Points, Polygon, Text, VLine};
use std::path::{Path, PathBuf};
use std::sync::{
    mpsc::{channel, Receiver, Sender},
    Arc,
};

const VERSION: &str = "0.4.0";
const WEIGHTS: [&str; 3] = ["none", "1/x", "1/x²"];
const README: &str = include_str!("../README.md");
const MAX_UNDO: usize = 10;

/// Scientific-notation styling for every plot (axes + hover text).
macro_rules! sciplot {
    ($p:expr) => {
        $p.x_axis_formatter(|m, _| sci(m.value))
            .y_axis_formatter(|m, _| sci(m.value))
            .label_formatter(|name, v| {
                let head = if name.is_empty() { String::new() } else { format!("{}\n", name) };
                format!("{}x = {}\ny = {}", head, sci(v.x), sci(v.y))
            })
    };
}

/// Chromatogram plots: retention-time axis with two decimals, other axes in scientific notation, no default cursor label.
macro_rules! rtplot {
    ($p:expr) => {
        $p.x_axis_formatter(|m, _| format!("{:.2}", m.value))
            .y_axis_formatter(|m, _| sci(m.value))
            .show_x(false)
            .show_y(false)
    };
}

/// Draws the traces, finds the point nearest to the cursor and annotates it with m/z, RT, scan number and intensity.
macro_rules! draw_series {
    ($pu:ident, $series:expr, $hover:ident, $scan_of:ident) => {{
        for s in $series.iter() {
            let d: Vec<[f64; 2]> = s.data.iter().map(|p| [p[0], p[1] * s.scale]).collect();
            let mut l = Line::new(PlotPoints::from(d)).name(s.name.clone()).color(pal(s.ci));
            if s.cur { l = l.width(2.5_f32); }
            $pu.line(l);
        }
        if let Some(p) = $pu.pointer_coordinate() {
            let b = $pu.plot_bounds();
            $hover = nearest_hover($series, p.x, p.y, b.max()[0] - b.min()[0]);
        }
        if let Some(h) = &$hover {
            let b = $pu.plot_bounds();
            let right = h.rt > (b.min()[0] + b.max()[0]) / 2.0;
            let txt = format!("{}\n{}\nRT {:.2} min · scan {}\nintensity {}", h.fname, h.desc, h.rt, $scan_of(h.file, h.level, h.idx), sci(h.y));
            $pu.points(Points::new(vec![[h.rt, h.y * h.scale]]).radius(5.0_f32).color(Color32::WHITE));
            $pu.text(Text::new(PlotPoint::new(h.rt, h.y * h.scale),
                RichText::new(txt).color(Color32::WHITE).background_color(Color32::from_black_alpha(200)))
                .anchor(if right { egui::Align2::RIGHT_BOTTOM } else { egui::Align2::LEFT_BOTTOM }));
        }
    }};
}

struct Series {
    name: String, data: Vec<[f64; 2]>, ci: usize, cur: bool, file: usize,
    /// 1 = points index the MS1 scan list, 2 = the MS2 scan list
    level: u8, fname: String, desc: String, scale: f64,
}

struct Hover { file: usize, rt: f64, y: f64, scale: f64, idx: usize, level: u8, fname: String, desc: String }

/// Point of any series closest to the cursor: nearest in time, then the trace whose intensity is closest to the cursor.
fn nearest_hover(series: &[Series], px: f64, py: f64, width: f64) -> Option<Hover> {
    let mut best: Option<(f64, Hover)> = None;
    for s in series {
        if s.data.is_empty() { continue; }
        let j = s.data.partition_point(|p| p[0] < px);
        let (c0, c1) = (j.saturating_sub(1), j.min(s.data.len() - 1));
        let bj = if (s.data[c1][0] - px).abs() < (s.data[c0][0] - px).abs() { c1 } else { c0 };
        let p = s.data[bj];
        if width > 0.0 && (p[0] - px).abs() > 0.05 * width { continue; }
        let dy = (p[1] * s.scale - py).abs();
        if best.as_ref().map(|b| dy < b.0).unwrap_or(true) {
            best = Some((dy, Hover { file: s.file, rt: p[0], y: p[1], scale: s.scale, idx: bj, level: s.level, fname: s.fname.clone(), desc: s.desc.clone() }));
        }
    }
    best.map(|b| b.1)
}

struct LoadInfo { path: String, raw_sha: String, detail: String }
enum Msg { Loaded(MzData, LoadInfo), Err(String), Finished }

const CREDITS: &str = "Thermo .raw files are read through ThermoRawFileParser (CompOmics, Apache-2.0), run as an external program; nothing is bundled or linked: https://github.com/compomics/ThermoRawFileParser\n\
Hulstaert N, Shofstahl J, Sachsenberg T, Walzer M, Barsnes H, Martens L, Perez-Riverol Y. J Proteome Res 2020;19:537-542 (doi:10.1021/acs.jproteome.9b00328).\n\
It uses Thermo Fisher Scientific's RawFileReader library (Copyright Thermo Fisher Scientific, Inc.), subject to Thermo's license.\n\
The external-converter approach follows mzmine: https://github.com/mzmine/mzmine (MIT).\n\
GUI: egui / eframe / egui_plot (MIT or Apache-2.0), https://github.com/emilk/egui; file dialogs: rfd (MIT).\n\
Built with Claude Sonnet 5.5 (Medium); app icon generated with Gemini.";

#[derive(Clone, Copy, PartialEq)]
enum Kind { Blank, Standard, QC, Unknown }
impl Kind {
    const ALL: [Kind; 4] = [Kind::Blank, Kind::Standard, Kind::QC, Kind::Unknown];
    fn label(self) -> &'static str {
        match self { Kind::Blank => "Blank", Kind::Standard => "Standard", Kind::QC => "QC", Kind::Unknown => "Unknown" }
    }
    fn from_label(s: &str) -> Kind {
        match s { "Blank" => Kind::Blank, "Standard" => Kind::Standard, "QC" => Kind::QC, _ => Kind::Unknown }
    }
}
#[derive(Clone, Copy, PartialEq)]
enum Mode { Areas, Curve }
#[derive(Clone, Copy, PartialEq)]
enum View { Xic, Tic, Bpc }
#[derive(Clone, Copy, PartialEq)]
enum Scope { Ions, AllTargets, Files }
#[derive(Clone, Copy, PartialEq)]
enum Tab { Chrom, Matrix }
#[derive(Clone, Copy, PartialEq)]
enum Basis { Ms2, Ms1 }

#[derive(Clone)]
struct Cell { start: f64, end: f64, found: bool, manual: bool, areas: Vec<f64>, logged: (f64, f64), include: bool }
impl Cell {
    fn new() -> Cell { Cell { start: 0.0, end: 0.0, found: false, manual: false, areas: vec![], logged: (0.0, 0.0), include: true } }
}

struct Sample {
    data: Arc<MzData>, src_path: String, raw_sha: String, kind: Kind, conc: f64,
    xics: Vec<Vec<Vec<[f64; 2]>>>, cells: Vec<Cell>,
    show: bool, rt_shift: f64, shift_ui: f64,
}

/// Everything undo/redo needs. Raw data is shared (Arc), derived chromatograms are recomputed on restore.
#[derive(Clone)]
struct SampleState {
    data: Arc<MzData>, src_path: String, raw_sha: String, kind: Kind, conc: f64,
    cells: Vec<Cell>, show: bool, rt_shift: f64,
}
#[derive(Clone)]
struct Snapshot {
    label: String, analytes: Vec<Analyte>, samples: Vec<SampleState>, sel_s: usize, sel_a: usize,
    params: Vec<(String, String)>,
}

struct Pending { recs: Vec<SampleRec>, mismatches: Vec<String> }

pub struct App {
    analytes: Vec<Analyte>,
    samples: Vec<Sample>,
    fits: Vec<Option<Fit>>,
    sel_s: usize,
    sel_a: usize,
    anchor: usize,
    tol_ppm: bool,
    tol: f64,
    rt_win: f64,
    weighting: u8,
    baseline: bool,
    ratio_tol: f64,
    units: String,
    reason: String,
    settings_dirty: bool,
    path_in: String,
    out_path: String,
    status: String,
    show_audit: bool,
    show_credits: bool,
    show_readme: bool,
    show_targets: bool,
    grab: Option<bool>,
    editing: bool,
    audit: Audit,
    user_logged: String,
    tx: Sender<Msg>,
    rx: Receiver<Msg>,
    loading: usize,
    trfp_path: String,
    mode: Mode,
    view: View,
    scope: Scope,
    tab: Tab,
    tic_level: u8,
    normalize: bool,
    matrix_conc: bool,
    align_on: bool,
    align_max: f64,
    align_ref: usize,
    history: History<Snapshot>,
    last_undo_label: String,
    rows: Vec<[String; 4]>,
    rows_dirty: bool,
    table_path: Option<PathBuf>,
    project_path: Option<PathBuf>,
    pending: Option<Pending>,
    zoom: f32,
    title: String,
    basis: Basis,
    show_spectrum: bool,
    spec_autofit: bool,
    follow_hover: bool,
    pinned_scan: Option<(usize, usize)>,
    hover: Option<(usize, f64)>,
    col_w: [f32; 4],
}

fn csv(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') { format!("\"{}\"", s.replace('"', "\"\"")) } else { s.to_string() }
}

fn calc_areas(xics: &[Vec<[f64; 2]>], c: &mut Cell, baseline: bool) {
    c.areas = xics.iter().map(|x| if c.found { integrate(x, c.start, c.end, baseline) } else { 0.0 }).collect();
}

fn pal(i: usize) -> Color32 {
    const P: [(u8, u8, u8); 10] = [(31, 119, 180), (214, 39, 40), (44, 160, 44), (255, 127, 14), (148, 103, 189),
                                   (140, 86, 75), (227, 119, 194), (127, 127, 127), (188, 189, 34), (23, 190, 207)];
    let c = P[i % 10];
    Color32::from_rgb(c.0, c.1, c.2)
}

fn chrom(d: &MzData, level: u8, bpc: bool) -> &Vec<[f64; 2]> {
    match (level, bpc) { (1, false) => &d.tic1, (1, true) => &d.bpc1, (_, false) => &d.tic2, (_, true) => &d.bpc2 }
}

fn align_trace(d: &MzData) -> &Vec<[f64; 2]> { if d.tic1.len() >= 5 { &d.tic1 } else { &d.tic2 } }

/// DragValue that shows and parses scientific notation.
fn dv(v: &mut f64) -> egui::DragValue<'_> {
    egui::DragValue::new(v).custom_formatter(|n, _| sci(n)).custom_parser(parse_sci)
}

/// DragValue for retention times / durations in minutes: two decimals ("7.20").
fn dv_rt(v: &mut f64) -> egui::DragValue<'_> {
    egui::DragValue::new(v).custom_formatter(|n, _| format!("{:.2}", n)).custom_parser(parse_sci)
}

fn parse_row(r: &[String; 4]) -> Option<Analyte> {
    let name = r[0].trim().to_string();
    if name.is_empty() { return None; }
    let na = |t: &str| t.is_empty() || t == "-" || t.eq_ignore_ascii_case("n/a") || t.eq_ignore_ascii_case("na");
    let t = r[1].trim();
    let rt = if na(t) { None } else { Some(t.parse::<f64>().ok()?) };
    let t = r[2].trim();
    let precursor = if na(t) { None } else { Some(t.parse::<f64>().ok()?) };
    let mut ions = vec![];
    for tok in r[3].split(|c: char| c == ';' || c == ',' || c.is_whitespace()).filter(|s| !s.is_empty()) { ions.push(tok.parse::<f64>().ok()?); }
    if ions.is_empty() { return None; }
    Some(Analyte { name, rt, ions, precursor })
}

impl App {
    pub fn new() -> App { App::new_with_audit("lcms_audit.log") }

    pub fn new_with_audit(audit_path: &str) -> App {
        let (tx, rx) = channel();
        let mut app = App {
            analytes: vec![], samples: vec![], fits: vec![], sel_s: 0, sel_a: 0, anchor: 0,
            tol_ppm: true, tol: 10.0, rt_win: 0.3, weighting: 2, baseline: true, ratio_tol: 30.0,
            units: "ng/mL".into(), reason: String::new(), settings_dirty: false,
            path_in: String::new(), out_path: "results.csv".into(),
            status: "Drop a .csv target table (optional) and Thermo .raw files onto the window, or use the File menu.".into(),
            show_audit: false, show_credits: false, show_readme: false, show_targets: false,
            grab: None, editing: false, audit: Audit::new(audit_path), user_logged: "N/A".into(), tx, rx, loading: 0,
            trfp_path: String::new(),
            mode: Mode::Curve, view: View::Xic, scope: Scope::Ions, tab: Tab::Chrom, tic_level: 1, normalize: false, matrix_conc: true,
            align_on: false, align_max: 1.0, align_ref: 0,
            history: History::new(MAX_UNDO), last_undo_label: String::new(),
            rows: vec![], rows_dirty: true, table_path: None, project_path: None, pending: None,
            zoom: 1.0, title: String::new(),
            basis: Basis::Ms2, show_spectrum: true, spec_autofit: true, follow_hover: true, pinned_scan: None, hover: None,
            col_w: [0.0; 4],
        };
        app.load_settings();
        app
    }

    // ---------------- settings (per machine) ----------------
    fn load_settings(&mut self) {
        if let Ok(s) = std::fs::read_to_string("lcms_settings.txt") {
            for l in s.lines() {
                if let Some(v) = l.strip_prefix("trfp=") { self.trfp_path = v.to_string(); }
                if let Some(v) = l.strip_prefix("user=") { if !v.trim().is_empty() { self.audit.user = v.trim().to_string(); self.user_logged = self.audit.user.clone(); } }
                if let Some(v) = l.strip_prefix("zoom=") { if let Ok(z) = v.parse::<f32>() { self.zoom = z.clamp(0.5, 3.0); } }
            }
        }
    }
    fn write_settings(&self) {
        let _ = std::fs::write("lcms_settings.txt", format!("trfp={}\nuser={}\nzoom={}\n", self.trfp_path.trim(), self.audit.user.trim(), self.zoom));
    }
    fn save_trfp(&mut self) {
        self.write_settings();
        self.audit.log("CONVERTER_PATHS", &format!("trfp={}", self.trfp_path.trim()));
    }

    fn tol_da(&self, mz: f64) -> f64 { if self.tol_ppm { mz * self.tol * 1e-6 } else { self.tol } }

    // ---------------- undo / redo ----------------
    fn params_vec(&self) -> Vec<(String, String)> {
        let kv = |k: &str, v: String| (k.to_string(), v);
        vec![
            kv("tol_ppm", self.tol_ppm.to_string()), kv("tol", self.tol.to_string()), kv("rt_win", self.rt_win.to_string()),
            kv("weighting", self.weighting.to_string()), kv("baseline", self.baseline.to_string()), kv("ratio_tol", self.ratio_tol.to_string()),
            kv("mode", if self.mode == Mode::Curve { "curve".into() } else { "areas".into() }), kv("units", self.units.clone()), kv("basis", if self.basis == Basis::Ms1 { "ms1".into() } else { "ms2".into() }),
            kv("align_on", self.align_on.to_string()), kv("align_max", self.align_max.to_string()), kv("align_ref", self.align_ref.to_string()),
        ]
    }
    fn apply_params(&mut self, kv: &[(String, String)]) {
        for (k, v) in kv {
            match k.as_str() {
                "tol_ppm" => if let Ok(x) = v.parse() { self.tol_ppm = x },
                "tol" => if let Ok(x) = v.parse() { self.tol = x },
                "rt_win" => if let Ok(x) = v.parse() { self.rt_win = x },
                "weighting" => if let Ok(x) = v.parse::<u8>() { self.weighting = x.min(2) },
                "baseline" => if let Ok(x) = v.parse() { self.baseline = x },
                "ratio_tol" => if let Ok(x) = v.parse() { self.ratio_tol = x },
                "mode" => self.mode = if v == "areas" { Mode::Areas } else { Mode::Curve },
                "units" => self.units = v.clone(),
                "basis" => self.basis = if v == "ms1" { Basis::Ms1 } else { Basis::Ms2 },
                "align_on" => if let Ok(x) = v.parse() { self.align_on = x },
                "align_max" => if let Ok(x) = v.parse() { self.align_max = x },
                "align_ref" => if let Ok(x) = v.parse() { self.align_ref = x },
                _ => {}
            }
        }
    }

    fn snapshot(&self, label: &str) -> Snapshot {
        Snapshot {
            label: label.to_string(), analytes: self.analytes.clone(),
            samples: self.samples.iter().map(|s| SampleState {
                data: s.data.clone(), src_path: s.src_path.clone(), raw_sha: s.raw_sha.clone(), kind: s.kind, conc: s.conc,
                cells: s.cells.clone(), show: s.show, rt_shift: s.rt_shift,
            }).collect(),
            sel_s: self.sel_s, sel_a: self.sel_a, params: self.params_vec(),
        }
    }

    fn restore(&mut self, s: Snapshot) {
        self.analytes = s.analytes;
        self.apply_params(&s.params);
        self.samples = s.samples.into_iter().map(|st| Sample {
            data: st.data, src_path: st.src_path, raw_sha: st.raw_sha, kind: st.kind, conc: st.conc, xics: vec![], cells: st.cells,
            show: st.show, rt_shift: st.rt_shift, shift_ui: st.rt_shift,
        }).collect();
        self.sel_s = s.sel_s.min(self.samples.len().saturating_sub(1));
        self.sel_a = s.sel_a.min(self.analytes.len().saturating_sub(1));
        self.anchor = self.sel_s;
        self.rows_dirty = true;
        self.editing = false;
        self.pinned_scan = None;
        self.recompute_all();
        self.fix_selection();
    }

    fn push_undo(&mut self, label: &str) {
        let s = self.snapshot(label);
        self.history.record(s);
        self.last_undo_label = label.to_string();
    }

    fn undo(&mut self) {
        let label = match self.history.peek_undo() { Some(s) => s.label.clone(), None => return };
        let cur = self.snapshot(&label);
        if let Some(prev) = self.history.undo(cur) {
            self.restore(prev);
            self.status = format!("Undid: {}", label);
            self.audit.log("UNDO", &label);
        }
    }

    fn redo(&mut self) {
        let label = match self.history.peek_redo() { Some(s) => s.label.clone(), None => return };
        let cur = self.snapshot(&label);
        if let Some(next) = self.history.redo(cur) {
            self.restore(next);
            self.status = format!("Redid: {}", label);
            self.audit.log("REDO", &label);
        }
    }

    // ---------------- selection ----------------
    fn cur(&self) -> Option<usize> {
        if self.sel_s < self.samples.len() && self.samples[self.sel_s].show { Some(self.sel_s) } else { None }
    }
    fn fix_selection(&mut self) {
        let n = self.samples.len();
        if n == 0 { self.sel_s = 0; self.anchor = 0; return; }
        if self.sel_s >= n { self.sel_s = n - 1; }
        if self.anchor >= n { self.anchor = n - 1; }
        if !self.samples[self.sel_s].show {
            if let Some(k) = (0..n).rev().find(|&k| self.samples[k].show) { self.sel_s = k; }
        }
    }
    fn shown(&self) -> Vec<usize> { (0..self.samples.len()).filter(|&i| self.samples[i].show).collect() }

    /// Explorer-style: click = only this, Ctrl+click = toggle, Shift+click = range from anchor (Ctrl+Shift adds the range).
    fn apply_click(&mut self, i: usize, ctrl: bool, shift: bool) {
        if i >= self.samples.len() { return; }
        if shift {
            let a = self.anchor.min(self.samples.len() - 1);
            let (lo, hi) = (a.min(i), a.max(i));
            if !ctrl { for s in self.samples.iter_mut() { s.show = false; } }
            for k in lo..=hi { self.samples[k].show = true; }
            self.sel_s = i;
        } else if ctrl {
            self.samples[i].show = !self.samples[i].show;
            self.anchor = i;
            if self.samples[i].show { self.sel_s = i; }
        } else {
            for s in self.samples.iter_mut() { s.show = false; }
            self.samples[i].show = true;
            self.sel_s = i;
            self.anchor = i;
        }
        self.fix_selection();
    }

    fn remove_sample(&mut self, i: usize) {
        if i >= self.samples.len() { return; }
        self.push_undo("remove file");
        let s = self.samples.remove(i);
        self.pinned_scan = None;
        self.audit.log("DATA_REMOVED", &format!("{} sha256={} (file on disk untouched)", s.data.name, s.data.sha256));
        if self.align_ref > i { self.align_ref -= 1; } else if self.align_ref == i { self.align_ref = 0; }
        self.fix_selection();
        for a in 0..self.analytes.len() { self.refit(a); }
        self.status = format!("Removed {} from the project (file on disk untouched)", s.data.name);
    }

    // ---------------- loading ----------------
    fn start_load(&mut self, ctx: &egui::Context, path: &str) -> bool {
        let raw = PathBuf::from(path);
        if !raw.is_file() { self.status = format!("Only single Thermo .raw files are supported / file not found: {}", path); return false; }
        let trfp = match rawconv::find_trfp(&self.trfp_path) {
            Some(p) => p,
            None => { self.status = "ThermoRawFileParser not found: enter the path to ThermoRawFileParser.exe in the top bar.".into(); return false; }
        };
        self.loading += 1;
        self.status = format!("Reading {} (held in memory only)...", path);
        let (tx, ctx, p) = (self.tx.clone(), ctx.clone(), path.to_string());
        std::thread::spawn(move || {
            let res = rawconv::hash_file(&raw).and_then(|h| rawconv::convert_in_memory(&trfp, &raw).map(|c| (h, c)));
            match res {
                Ok((h, c)) => {
                    let name = raw.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
                    match mzml::parse_bytes(&name, &c.mzml) {
                        Ok(d) => { let _ = tx.send(Msg::Loaded(d, LoadInfo { path: p, raw_sha: h, detail: c.detail })); }
                        Err(e) => { let _ = tx.send(Msg::Err(e)); }
                    }
                }
                Err(e) => { let _ = tx.send(Msg::Err(format!("{}: {}", p, e))); }
            }
            let _ = tx.send(Msg::Finished);
            ctx.request_repaint();
        });
        true
    }

    fn load_path(&mut self, ctx: &egui::Context, path: &str) {
        let lower = path.to_lowercase();
        if lower.ends_with(".raw") {
            self.start_load(ctx, path);
        } else if lower.ends_with(".quant") {
            self.open_project(ctx, Path::new(path));
        } else if lower.ends_with(".csv") || lower.ends_with(".tsv") || lower.ends_with(".txt") {
            self.load_table(Path::new(path));
        } else {
            self.status = "Unsupported file: drop Thermo .raw files, a .csv target table or a .quant project.".into();
        }
    }

    fn load_table(&mut self, path: &Path) {
        match std::fs::read(path) {
            Ok(bytes) => match parse_table(&String::from_utf8_lossy(&bytes)) {
                Ok(a) => {
                    self.push_undo("load target table");
                    self.audit.log("TABLE_LOADED", &format!("{} sha256={} analytes={}", path.display(), crate::sha256::sha256_hex(&bytes), a.len()));
                    self.analytes = a;
                    self.sel_a = 0;
                    self.table_path = Some(path.to_path_buf());
                    self.rows_dirty = true;
                    for s in self.samples.iter_mut() { s.cells.clear(); }
                    self.recompute_all();
                    self.status = format!("Loaded {} targets", self.analytes.len());
                }
                Err(e) => self.status = format!("Table error: {}", e),
            },
            Err(e) => self.status = format!("{}: {}", path.display(), e),
        }
    }

    fn add_sample(&mut self, d: MzData, info: LoadInfo) {
        self.audit.log("DATA_LOADED", &format!("{} sha256={} ms2_scans={} ms1_scans={} total_scans={} source_raw={} raw_sha256={} {}",
            d.name, d.sha256, d.spectra.len(), d.tic1.len(), d.n_total, info.path, info.raw_sha, info.detail));
        let in_project = self.pending.is_some();
        if !in_project { self.push_undo("import file"); }
        let na = self.analytes.len();
        let mut sm = Sample {
            data: Arc::new(d), src_path: info.path.clone(), raw_sha: info.raw_sha.clone(), kind: Kind::Unknown, conc: 0.0,
            xics: vec![], cells: vec![], show: true, rt_shift: 0.0, shift_ui: 0.0,
        };
        if let Some(p) = self.pending.as_mut() {
            if let Some(rec) = p.recs.iter().find(|r| r.path == info.path) {
                sm.kind = Kind::from_label(&rec.kind);
                sm.conc = rec.conc;
                sm.rt_shift = rec.rt_shift;
                sm.shift_ui = rec.rt_shift;
                sm.show = rec.show;
                sm.cells = (0..na).map(|_| Cell::new()).collect();
                for c in &rec.cells {
                    if c.analyte < na {
                        let cell = &mut sm.cells[c.analyte];
                        cell.found = c.found; cell.manual = c.manual; cell.start = c.start; cell.end = c.end;
                        cell.logged = (c.start, c.end); cell.include = c.include;
                    }
                }
                if !rec.raw_sha256.is_empty() && rec.raw_sha256 != info.raw_sha { p.mismatches.push(sm.data.name.clone()); }
            }
        }
        self.samples.push(sm);
        let i = self.samples.len() - 1;
        if !in_project {
            self.sel_s = i;
            self.anchor = i;
            if self.align_on && self.align_ref < i {
                if let Some(sh) = self.compute_shift(i) {
                    self.samples[i].rt_shift = sh;
                    self.samples[i].shift_ui = sh;
                    let msg = format!("file={} sha256={} reference={} shift_min={:.4} max_search_min={}",
                        self.samples[i].data.name, self.samples[i].data.sha256, self.samples[self.align_ref].data.name, sh, self.align_max);
                    self.audit.log("RT_ALIGNMENT", &msg);
                }
            }
        }
        self.recompute_sample(i);
        for a in 0..self.analytes.len() { self.refit(a); }
        self.status = format!("Loaded {}", self.samples[i].data.name);
    }

    // ---------------- computation ----------------
    fn recompute_sample(&mut self, s: usize) {
        let (na, rt_win, baseline) = (self.analytes.len(), self.rt_win, self.baseline);
        let basis = self.basis;
        let traces: Vec<Vec<f64>> = (0..na).map(|a| self.trace_mz(a)).collect();
        let tols: Vec<Vec<f64>> = traces.iter().map(|v| v.iter().map(|&m| self.tol_da(m)).collect()).collect();
        let sm = &mut self.samples[s];
        let shift = sm.rt_shift;
        sm.xics = traces.iter().enumerate()
            .map(|(ai, mzs)| mzs.iter().enumerate().map(|(ii, &mz)| {
                let sp: &[Spectrum] = if basis == Basis::Ms1 { &sm.data.ms1 } else { &sm.data.spectra };
                let mut x = xic(sp, mz, tols[ai][ii]);
                if shift != 0.0 { for p in x.iter_mut() { p[0] += shift; } }
                x
            }).collect()).collect();
        if sm.cells.len() != na { sm.cells = (0..na).map(|_| Cell::new()).collect(); }
        for a in 0..na {
            let c = &mut sm.cells[a];
            if !c.manual {
                match auto_peak(&sm.xics[a][0], self.analytes[a].rt, rt_win) {
                    Some((st, en)) => { c.start = st; c.end = en; c.found = true; c.logged = (st, en); }
                    None => c.found = false,
                }
            }
            calc_areas(&sm.xics[a], c, baseline);
        }
    }

    /// Recompute a single (new or edited) target in every file; its integrations are reset.
    fn recompute_analyte(&mut self, a: usize) {
        let na = self.analytes.len();
        if a >= na { return; }
        let (rt_win, baseline) = (self.rt_win, self.baseline);
        let basis = self.basis;
        let mzs = self.trace_mz(a);
        let tols: Vec<f64> = mzs.iter().map(|&m| self.tol_da(m)).collect();
        for s in 0..self.samples.len() {
            let sm = &mut self.samples[s];
            if sm.xics.len() < na { sm.xics.resize(na, vec![]); }
            if sm.cells.len() < na { sm.cells.resize_with(na, Cell::new); }
            let shift = sm.rt_shift;
            sm.xics[a] = mzs.iter().enumerate().map(|(ii, &mz)| {
                let sp: &[Spectrum] = if basis == Basis::Ms1 { &sm.data.ms1 } else { &sm.data.spectra };
                let mut x = xic(sp, mz, tols[ii]);
                if shift != 0.0 { for p in x.iter_mut() { p[0] += shift; } }
                x
            }).collect();
            let mut c = Cell::new();
            if let Some((st, en)) = auto_peak(&sm.xics[a][0], self.analytes[a].rt, rt_win) {
                c.start = st; c.end = en; c.found = true; c.logged = (st, en);
            }
            calc_areas(&sm.xics[a], &mut c, baseline);
            sm.cells[a] = c;
        }
        if self.fits.len() < na { self.fits.resize(na, None); }
        self.refit(a);
    }

    fn recompute_all(&mut self) {
        self.fits = vec![None; self.analytes.len()];
        for s in 0..self.samples.len() { self.recompute_sample(s); }
        for a in 0..self.analytes.len() { self.refit(a); }
    }

    fn refit(&mut self, a: usize) {
        if a >= self.analytes.len() { return; }
        if self.fits.len() != self.analytes.len() { self.fits = vec![None; self.analytes.len()]; }
        let pts: Vec<(f64, f64)> = self.samples.iter()
            .filter(|s| s.kind == Kind::Standard && s.conc > 0.0 && s.cells.len() > a && s.cells[a].found && s.cells[a].include)
            .map(|s| (s.conc, s.cells[a].areas[0])).collect();
        self.fits[a] = fit(&pts, self.weighting);
    }

    fn apply_edit(&mut self, s: usize, a: usize, st: f64, en: f64) {
        if !self.editing { self.push_undo("manual integration"); self.editing = true; }
        let baseline = self.baseline;
        let sm = &mut self.samples[s];
        let c = &mut sm.cells[a];
        c.start = st; c.end = en; c.found = true; c.manual = true;
        calc_areas(&sm.xics[a], c, baseline);
        self.refit(a);
    }

    fn commit(&mut self, s: usize, a: usize) {
        self.editing = false;
        let (st, en, area, logged) = {
            let c = &self.samples[s].cells[a];
            (c.start, c.end, c.areas.first().copied().unwrap_or(0.0), c.logged)
        };
        if logged != (st, en) {
            let msg = format!("file={} sha256={} analyte={} start {:.4}->{:.4} end {:.4}->{:.4} area={:.1} reason={}",
                self.samples[s].data.name, self.samples[s].data.sha256, self.analytes[a].name,
                logged.0, st, logged.1, en, area,
                if self.reason.trim().is_empty() { "NONE GIVEN" } else { self.reason.trim() });
            self.samples[s].cells[a].logged = (st, en);
            self.audit.log("MANUAL_INTEGRATION", &msg);
        }
    }

    fn apply_settings(&mut self) {
        self.audit.log("PARAMETERS", &format!("tol={}{} rt_window={} weighting={} baseline={} ion_ratio_tol={}%",
            self.tol, if self.tol_ppm { "ppm" } else { "Da" }, self.rt_win, WEIGHTS[self.weighting as usize], self.baseline, self.ratio_tol));
        self.recompute_all();
        self.settings_dirty = false;
    }

    // ---------------- retention-time alignment ----------------
    fn compute_shift(&self, s: usize) -> Option<f64> {
        let r = self.samples.get(self.align_ref)?;
        let x = &self.samples[s];
        let (tr, tx) = (align_trace(&r.data), align_trace(&x.data));
        if tr.len() < 5 || tx.len() < 5 { return None; }
        Some(estimate_shift(tr, tx, self.align_max) + r.rt_shift)
    }

    fn set_shift(&mut self, s: usize, shift: f64, why: &str) {
        let old = self.samples[s].rt_shift;
        {
            let sm = &mut self.samples[s];
            sm.rt_shift = shift;
            sm.shift_ui = shift;
            for c in sm.cells.iter_mut() { c.manual = false; }
        }
        self.recompute_sample(s);
        for a in 0..self.analytes.len() { self.refit(a); }
        let msg = format!("file={} sha256={} shift {:.4}->{:.4} min reason={} (manual integrations of this file reset)",
            self.samples[s].data.name, self.samples[s].data.sha256, old, shift, why);
        self.audit.log("RT_SHIFT", &msg);
    }

    fn align_all(&mut self) {
        let n = self.samples.len();
        if n < 2 { self.status = "Need at least two files to align.".into(); return; }
        if self.align_ref >= n { self.align_ref = 0; }
        self.push_undo("align all");
        for s in 0..n {
            if s == self.align_ref { continue; }
            if let Some(sh) = self.compute_shift(s) { self.set_shift(s, sh, "align all (TIC cross-correlation)"); }
        }
        self.status = format!("Aligned {} files to {}", n - 1, self.samples[self.align_ref].data.name);
    }

    fn clear_alignment(&mut self) {
        self.push_undo("clear alignment");
        for s in 0..self.samples.len() { if self.samples[s].rt_shift != 0.0 { self.set_shift(s, 0.0, "clear alignment"); } }
        self.status = "Alignment cleared".into();
    }

    // ---------------- results ----------------
    fn conc_of(&self, s: usize, a: usize) -> Option<f64> {
        if self.mode == Mode::Areas { return None; }
        let f = self.fits.get(a).copied().flatten()?;
        let c = self.samples.get(s)?.cells.get(a)?;
        if !c.found || f.m == 0.0 { return None; }
        Some((c.areas[0] - f.b) / f.m)
    }

    fn ref_ratio(&self, a: usize, ion: usize) -> Option<f64> {
        let v: Vec<f64> = self.samples.iter()
            .filter(|s| s.kind == Kind::Standard && s.cells.len() > a && s.cells[a].found && s.cells[a].include && s.cells[a].areas[0] > 0.0)
            .map(|s| s.cells[a].areas[ion] / s.cells[a].areas[0]).collect();
        if v.is_empty() { None } else { Some(v.iter().sum::<f64>() / v.len() as f64) }
    }

    fn flags(&self, s: usize, a: usize) -> String {
        let sm = &self.samples[s];
        let c = match sm.cells.get(a) { Some(c) => c, None => return String::new() };
        let mut f: Vec<String> = vec![];
        if !c.found || c.areas[0] <= 0.0 { return "ND".into(); }
        if c.manual { f.push("manual".into()); }
        if sm.kind == Kind::Blank { f.push("PEAK IN BLANK".into()); }
        if let (Some(fit), Some(cc)) = (self.fits.get(a).copied().flatten(), self.conc_of(s, a)) {
            if sm.kind != Kind::Standard && (cc < fit.lo || cc > fit.hi) { f.push("outside cal range".into()); }
            if sm.kind == Kind::Standard && sm.conc > 0.0 && ((cc - sm.conc) / sm.conc).abs() > 0.15 { f.push("accuracy >15%".into()); }
        }
        for i in 1..c.areas.len() {
            if let Some(r) = self.ref_ratio(a, i) {
                let obs = c.areas[i] / c.areas[0];
                if r > 0.0 && ((obs - r) / r).abs() * 100.0 > self.ratio_tol { f.push(format!("ion{} ratio fail", i + 1)); }
            }
        }
        f.join("; ")
    }

    fn export(&mut self) -> Result<(), String> {
        let curve = self.mode == Mode::Curve;
        let mut o = String::new();
        o += &format!("# LC-MS/MS Quant v{} export {} user={}\n", VERSION, iso_now(), self.audit.user);
        o += &format!("# quantification={}\n", if curve { "standard curve" } else { "integrated areas only (arbitrary units)" });
        o += &format!("# quantified_signal={}\n", if self.basis == Basis::Ms1 { "MS1 precursor area" } else { "MS2 product-ion area" });
        o += &format!("# tol={}{} rt_window_min={} weighting={} baseline_subtraction={} ion_ratio_tol_pct={} conc_units={}\n",
            self.tol, if self.tol_ppm { "ppm" } else { "Da" }, self.rt_win, WEIGHTS[self.weighting as usize], self.baseline, self.ratio_tol, self.units);
        o += &format!("# rt_alignment_on_import={} max_shift_min={} reference={}\n", self.align_on, self.align_max,
            self.samples.get(self.align_ref).map(|s| s.data.name.clone()).unwrap_or_default());
        if curve {
            for (a, an) in self.analytes.iter().enumerate() {
                match self.fits.get(a).copied().flatten() {
                    Some(f) => o += &format!("# CAL {}: area = {:.6e}*conc + {:.6e}; R2={:.5}; n={}; range {}-{}\n", an.name, f.m, f.b, f.r2, f.n, f.lo, f.hi),
                    None => o += &format!("# CAL {}: no valid calibration\n", an.name),
                }
            }
        }
        o += "file,sha256,type,level,rt_shift_min,analyte,start_min,end_min,area_quant_(int*s),manual,ion_ratios,calc_conc,flags\n";
        for s in 0..self.samples.len() {
            for a in 0..self.analytes.len() {
                let sm = &self.samples[s];
                let c = &sm.cells[a];
                let ratios: Vec<String> = (1..c.areas.len()).map(|i| if c.areas[0] > 0.0 { format!("{:.3}", c.areas[i] / c.areas[0]) } else { "n/a".into() }).collect();
                let conc = self.conc_of(s, a).map(|v| format!("{:.4}", v)).unwrap_or("n/a".into());
                o += &format!("{},{},{},{},{:.4},{},{:.4},{:.4},{:.1},{},{},{},{}\n", csv(&sm.data.name), sm.data.sha256, sm.kind.label(),
                    if sm.kind == Kind::Standard && curve { sm.conc.to_string() } else { String::new() }, sm.rt_shift, csv(&self.analytes[a].name),
                    c.start, c.end, c.areas[0], c.manual, ratios.join(" "), conc, csv(&self.flags(s, a)));
            }
        }
        std::fs::write(&self.out_path, o).map_err(|e| e.to_string())?;
        self.audit.log("EXPORT", &self.out_path.clone());
        Ok(())
    }

    fn choose_out_path(&mut self) {
        let mut d = rfd::FileDialog::new().set_title("Select analysis output file").add_filter("CSV", &["csv"]);
        let p = PathBuf::from(&self.out_path);
        if let Some(dir) = p.parent() { if dir.is_dir() { d = d.set_directory(dir); } }
        d = d.set_file_name(p.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_else(|| "results.csv".into()));
        if let Some(path) = d.save_file() { self.out_path = path.to_string_lossy().to_string(); }
    }

    // ---------------- target table editing ----------------
    fn sync_rows(&mut self) {
        self.rows = self.analytes.iter().map(|a| [
            a.name.clone(),
            a.rt.map(|v| v.to_string()).unwrap_or_else(|| "n/a".into()),
            a.precursor.map(|v| v.to_string()).unwrap_or_else(|| "n/a".into()),
            a.ions.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("; "),
        ]).collect();
        self.rows_dirty = false;
    }

    fn apply_row(&mut self, i: usize) {
        if i >= self.rows.len() || i >= self.analytes.len() { return; }
        if let Some(an) = parse_row(&self.rows[i]) {
            if an != self.analytes[i] {
                self.push_undo("edit target");
                self.audit.log("TARGET_EDIT", &format!("#{} {} -> {} rt={:?} ions={:?}", i + 1, self.analytes[i].name, an.name, an.rt, an.ions));
                self.analytes[i] = an;
                self.recompute_analyte(i);
            }
        }
    }

    fn add_target(&mut self) {
        self.push_undo("add target");
        let n = self.analytes.len();
        self.analytes.push(Analyte { name: format!("Target {}", n + 1), rt: None, ions: vec![100.0], precursor: None });
        self.fits.push(None);
        self.recompute_analyte(n);
        self.rows_dirty = true;
        self.audit.log("TARGET_ADDED", &format!("#{}", n + 1));
    }

    fn delete_target(&mut self, i: usize) {
        if i >= self.analytes.len() { return; }
        self.push_undo("delete target");
        let name = self.analytes.remove(i).name;
        for sm in self.samples.iter_mut() {
            if i < sm.xics.len() { sm.xics.remove(i); }
            if i < sm.cells.len() { sm.cells.remove(i); }
        }
        if i < self.fits.len() { self.fits.remove(i); }
        self.sel_a = self.sel_a.min(self.analytes.len().saturating_sub(1));
        self.rows_dirty = true;
        self.audit.log("TARGET_DELETED", &name);
    }

    fn targets_csv(&self) -> String {
        let mut o = String::from("name,rt_min,ions\n");
        for a in &self.analytes {
            let mut ions = a.ions.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(";");
            if let Some(p) = a.precursor { ions += &format!(";p:{}", p); }
            o += &format!("{},{},{}\n", csv(&a.name), a.rt.map(|v| v.to_string()).unwrap_or_else(|| "n/a".into()), ions);
        }
        o
    }

    fn write_targets(&mut self, path: &Path) {
        let t = self.targets_csv();
        match std::fs::write(path, &t) {
            Ok(_) => {
                self.audit.log("TARGETS_SAVED", &format!("{} sha256={}", path.display(), crate::sha256::sha256_hex(t.as_bytes())));
                self.table_path = Some(path.to_path_buf());
                self.status = format!("Saved target table {}", path.display());
            }
            Err(e) => self.status = format!("Save failed: {}", e),
        }
    }

    fn save_targets_as(&mut self) {
        let mut d = rfd::FileDialog::new().set_title("Save target table").add_filter("CSV", &["csv"]).set_file_name("targets.csv");
        if let Some(p) = &self.table_path {
            if let Some(dir) = p.parent() { d = d.set_directory(dir); }
            if let Some(n) = p.file_name() { d = d.set_file_name(n.to_string_lossy().to_string()); }
        }
        if let Some(path) = d.save_file() { self.write_targets(&path); }
    }

    fn open_targets_dialog(&mut self) {
        if let Some(p) = rfd::FileDialog::new().set_title("Open target table").add_filter("Target table", &["csv", "tsv", "txt"]).pick_file() { self.load_table(&p); }
    }

    // ---------------- projects ----------------
    fn project_file(&self) -> ProjectFile {
        ProjectFile {
            app_version: VERSION.into(), user: self.audit.user.clone(), params: self.params_vec(), analytes: self.analytes.clone(),
            samples: self.samples.iter().map(|s| SampleRec {
                path: s.src_path.clone(), raw_sha256: s.raw_sha.clone(), kind: s.kind.label().into(), conc: s.conc, rt_shift: s.rt_shift, show: s.show,
                cells: s.cells.iter().enumerate().map(|(i, c)| CellRec { analyte: i, found: c.found, manual: c.manual, start: c.start, end: c.end, include: c.include }).collect(),
            }).collect(),
        }
    }

    fn save_project(&mut self, force_dialog: bool) {
        let path = match (&self.project_path, force_dialog) {
            (Some(p), false) => p.clone(),
            _ => match rfd::FileDialog::new().set_title("Save project").add_filter("LC-MS/MS Quant project", &["quant"]).set_file_name("project.quant").save_file() {
                Some(p) => if p.extension().is_none() { p.with_extension("quant") } else { p },
                None => return,
            },
        };
        let text = self.project_file().to_text();
        match std::fs::write(&path, &text) {
            Ok(_) => {
                self.audit.log("PROJECT_SAVED", &format!("{} sha256={}", path.display(), crate::sha256::sha256_hex(text.as_bytes())));
                self.status = format!("Saved project {}", path.display());
                self.project_path = Some(path);
            }
            Err(e) => self.status = format!("Save failed: {}", e),
        }
    }

    fn open_project_dialog(&mut self, ctx: &egui::Context) {
        if let Some(p) = rfd::FileDialog::new().set_title("Open project").add_filter("LC-MS/MS Quant project", &["quant"]).pick_file() { self.open_project(ctx, &p); }
    }

    fn open_project(&mut self, ctx: &egui::Context, path: &Path) {
        let text = match std::fs::read_to_string(path) { Ok(t) => t, Err(e) => { self.status = format!("{}: {}", path.display(), e); return; } };
        let pf = match ProjectFile::parse(&text) { Ok(p) => p, Err(e) => { self.status = format!("Project error: {}", e); return; } };
        self.push_undo("open project");
        self.audit.log("PROJECT_OPENED", &format!("{} sha256={} saved_by={}", path.display(), crate::sha256::sha256_hex(text.as_bytes()), pf.user));
        self.samples.clear();
        self.analytes = pf.analytes.clone();
        self.fits = vec![None; self.analytes.len()];
        self.apply_params(&pf.params);
        self.sel_s = 0; self.sel_a = 0; self.anchor = 0;
        self.rows_dirty = true;
        self.project_path = Some(path.to_path_buf());
        if pf.samples.is_empty() { self.status = format!("Opened project {}", path.display()); return; }
        self.pending = Some(Pending { recs: pf.samples.clone(), mismatches: vec![] });
        let mut started = 0;
        for rec in &pf.samples { if self.start_load(ctx, &rec.path) { started += 1; } }
        if started == 0 { self.pending = None; }
        else { self.status = format!("Opening project: re-importing {} raw file(s)...", started); }
    }

    fn finish_project(&mut self, p: Pending) {
        let recs = &p.recs;
        self.samples.sort_by_key(|s| recs.iter().position(|r| r.path == s.src_path).unwrap_or(usize::MAX));
        let n = self.samples.len();
        self.fix_selection();
        for a in 0..self.analytes.len() { self.refit(a); }
        if self.align_ref >= n { self.align_ref = 0; }
        let mut msg = format!("Project loaded: {} of {} raw files", n, recs.len());
        if !p.mismatches.is_empty() {
            msg += &format!("  ⚠ SHA-256 differs from the saved project for: {}", p.mismatches.join(", "));
            self.audit.log("PROJECT_HASH_MISMATCH", &p.mismatches.join(", "));
        }
        if n < recs.len() { msg += "  ⚠ some files could not be re-imported (moved or converter missing)"; }
        self.status = msg;
    }

    fn new_project(&mut self) {
        self.push_undo("new project");
        self.audit.log("PROJECT_NEW", "");
        self.samples.clear();
        self.analytes.clear();
        self.fits.clear();
        self.sel_s = 0; self.sel_a = 0; self.anchor = 0;
        self.project_path = None;
        self.table_path = None;
        self.pending = None;
        self.rows_dirty = true;
        self.status = "New project".into();
    }

    // ---------------- UI ----------------
    /// Selected target, all traces, in file `s`, with draggable integration boundaries and hover read-out.
    fn ui_plot(&mut self, ui: &mut egui::Ui, s: usize) {
        let a = self.sel_a;
        let mut grab = self.grab;
        let mut edit: Option<(f64, f64)> = None;
        let mut stopped = false;
        let clicked;
        let mut hover: Option<Hover> = None;
        let tmz = self.trace_mz(a);
        {
            let sm = &self.samples[s];
            let c = &sm.cells[a];
            let norm = self.normalize;
            let level: u8 = if self.basis == Basis::Ms1 { 1 } else { 2 };
            let series: Vec<Series> = sm.xics[a].iter().enumerate().map(|(i, x)| {
                let mx = x.iter().fold(0.0f64, |m, p| m.max(p[1]));
                let d = format!("m/z {}", sci_prec(tmz[i], 9));
                Series { name: d.clone(), data: x.clone(), ci: i, cur: false, file: s, level, fname: sm.data.name.clone(),
                         desc: format!("{} · {}", self.analytes[a].name, d), scale: if norm && mx > 0.0 { 1.0 / mx } else { 1.0 } }
            }).collect();
            let scan_of = |f: usize, l: u8, j: usize| -> u32 {
                let d = &self.samples[f].data;
                if l == 1 { d.ms1.get(j).map(|x| x.scan).unwrap_or(0) } else { d.spectra.get(j).map(|x| x.scan).unwrap_or(0) }
            };
            let pr = rtplot!(Plot::new("xic").height(300.0).allow_drag(false).allow_boxed_zoom(false)
                .legend(Legend::default()).x_axis_label("Retention time (min)  |  drag a green/red boundary line, click = pin the scan"))
                .show(ui, |pu| {
                    draw_series!(pu, &series, hover, scan_of);
                    if c.found && !norm {
                        pu.polygon(Polygon::new(PlotPoints::from(shade(&sm.xics[a][0], c.start, c.end)))
                            .fill_color(Color32::from_rgba_unmultiplied(31, 119, 180, 60)).name("integrated"));
                    }
                    if c.found {
                        pu.vline(VLine::new(c.start).color(Color32::DARK_GREEN));
                        pu.vline(VLine::new(c.end).color(Color32::DARK_RED));
                    }
                    let r = pu.response().clone();
                    if let Some(p) = pu.pointer_coordinate() {
                        if r.drag_started_by(PointerButton::Primary) && c.found {
                            grab = Some((p.x - c.start).abs() <= (p.x - c.end).abs());
                        }
                        if r.dragged_by(PointerButton::Primary) {
                            if let Some(g) = grab {
                                let (mut st, mut en) = (c.start, c.end);
                                if g { st = p.x.min(en - 1e-4); } else { en = p.x.max(st + 1e-4); }
                                edit = Some((st, en));
                            }
                        }
                    }
                    if r.drag_stopped() { stopped = true; }
                });
            clicked = pr.response.clicked();
        }
        if let Some(h) = &hover {
            self.hover = Some((h.file, h.rt));
            if clicked { self.pin_scan(h.file, h.rt); }
        }
        if let Some((st, en)) = edit { self.apply_edit(s, a, st, en); }
        if stopped { self.commit(s, a); grab = None; }
        self.grab = grab;

        let (found, mut st, mut en) = { let c = &self.samples[s].cells[a]; (c.found, c.start, c.end) };
        let mut changed = false; let mut done = false;
        ui.horizontal(|ui| {
            if found {
                ui.label("Start (min)");
                let r1 = ui.add(dv_rt(&mut st).speed(0.01));
                ui.label("End (min)");
                let r2 = ui.add(dv_rt(&mut en).speed(0.01));
                changed = r1.changed() || r2.changed();
                done = r1.drag_stopped() || r1.lost_focus() || r2.drag_stopped() || r2.lost_focus();
            } else if ui.button("No peak found - define manually").clicked() {
                let x = &self.samples[s].xics[a][0];
                let mid = self.analytes[a].rt.unwrap_or_else(|| x.get(x.len() / 2).map(|p| p[0]).unwrap_or(1.0));
                st = mid - 0.1; en = mid + 0.1; changed = true; done = true;
            }
            ui.label("Reason for manual change:");
            ui.add(egui::TextEdit::singleline(&mut self.reason).desired_width(260.0));
        });
        if changed && st < en { self.apply_edit(s, a, st, en); }
        if done { self.commit(s, a); }
        if found && ui.button("Reset to automatic boundaries").clicked() {
            self.push_undo("reset integration");
            self.audit.log("MANUAL_RESET", &format!("file={} analyte={}", self.samples[s].data.name, self.analytes[a].name));
            self.samples[s].cells[a].manual = false;
            self.recompute_sample(s);
            self.refit(a);
        }
        let c = &self.samples[s].cells[a];
        if c.found {
            let mut t = String::new();
            for (i, ar) in c.areas.iter().enumerate() {
                t += &format!("{} {} (m/z {}): area {}   ", if self.basis == Basis::Ms1 { "precursor" } else { "ion" }, i + 1, sci_prec(tmz[i], 9), sci(*ar));
                if i > 0 && c.areas[0] > 0.0 { t += &format!("[ratio {}]   ", sci(ar / c.areas[0])); }
            }
            ui.monospace(t);
            match self.conc_of(s, a) {
                Some(v) => ui.strong(format!("Calculated: {} {}   {}", sci(v), self.units, self.flags(s, a))),
                None => if self.mode == Mode::Curve {
                    ui.label("No calibration available for this analyte")
                } else {
                    ui.strong(format!("Integrated area: {} a.u. (intensity·s)   {}", sci(c.areas[0]), self.flags(s, a)))
                },
            };
        }
    }

    /// Read-only overlays: all targets of the current file, one target across selected files, or TIC/BPC of selected files.
    fn ui_overlay(&mut self, ui: &mut egui::Ui, cur: Option<usize>) {
        let a = self.sel_a;
        let norm = self.normalize;
        let shown = self.shown();
        let level: u8 = if self.basis == Basis::Ms1 { 1 } else { 2 };
        let mut series: Vec<Series> = vec![];
        match (self.view, self.scope) {
            (View::Xic, Scope::AllTargets) => {
                if let Some(s) = cur {
                    let sm = &self.samples[s];
                    for (ai, an) in self.analytes.iter().enumerate() {
                        if let Some(x) = sm.xics.get(ai).and_then(|v| v.first()) {
                            let desc = format!("{} · m/z {}", an.name, sci_prec(self.trace_mz(ai)[0], 9));
                            series.push(Series { name: an.name.clone(), data: x.clone(), ci: ai, cur: false, file: s, level, fname: sm.data.name.clone(), desc, scale: 1.0 });
                        }
                    }
                }
            }
            (View::Xic, _) => {
                let desc = format!("{} · m/z {}", self.analytes[a].name, sci_prec(self.trace_mz(a)[0], 9));
                for &f in &shown {
                    let sm = &self.samples[f];
                    if let Some(x) = sm.xics.get(a).and_then(|v| v.first()) {
                        series.push(Series { name: sm.data.name.clone(), data: x.clone(), ci: f, cur: Some(f) == cur, file: f, level, fname: sm.data.name.clone(), desc: desc.clone(), scale: 1.0 });
                    }
                }
            }
            (v, _) => {
                let bpc = v == View::Bpc;
                for &f in &shown {
                    let sm = &self.samples[f];
                    let sh = sm.rt_shift;
                    let d: Vec<[f64; 2]> = chrom(&sm.data, self.tic_level, bpc).iter().map(|p| [p[0] + sh, p[1]]).collect();
                    series.push(Series { name: sm.data.name.clone(), data: d, ci: f, cur: Some(f) == cur, file: f, level: self.tic_level, fname: sm.data.name.clone(),
                                         desc: format!("{} (MS{})", if bpc { "BPC" } else { "TIC" }, self.tic_level), scale: 1.0 });
                }
            }
        }
        series.retain(|x| !x.data.is_empty());
        if shown.is_empty() { ui.label("No file selected."); return; }
        if series.is_empty() {
            ui.label(if self.view == View::Xic { "No chromatogram data." } else { "The selected file(s) have no scans at this MS level." });
            return;
        }
        if norm {
            for sr in series.iter_mut() {
                let m = sr.data.iter().fold(0.0f64, |m, p| m.max(p[1]));
                if m > 0.0 { sr.scale = 1.0 / m; }
            }
        }
        let rt_hint = if self.view == View::Xic && !self.analytes.is_empty() { self.analytes[a].rt } else { None };
        let id = (self.view as u8, self.scope as u8, self.tic_level);
        let mut hover: Option<Hover> = None;
        let clicked;
        {
            let scan_of = |f: usize, l: u8, j: usize| -> u32 {
                let d = &self.samples[f].data;
                if l == 1 { d.ms1.get(j).map(|x| x.scan).unwrap_or(0) } else { d.spectra.get(j).map(|x| x.scan).unwrap_or(0) }
            };
            let pr = rtplot!(Plot::new(id).height(380.0).legend(Legend::default()).x_axis_label("Retention time (min)  |  click = pin the nearest MS2 scan")
                .y_axis_label(if norm { "Relative intensity" } else { "Intensity" }))
                .show(ui, |pu| {
                    draw_series!(pu, &series, hover, scan_of);
                    if let Some(h) = rt_hint { pu.vline(VLine::new(h).color(Color32::GRAY).name("expected RT")); }
                });
            clicked = pr.response.clicked();
        }
        if let Some(h) = &hover {
            self.hover = Some((h.file, h.rt));
            if clicked { self.pin_scan(h.file, h.rt); }
        }
    }

    fn ui_shift(&mut self, ui: &mut egui::Ui, s: usize) {
        let mut v = self.samples[s].shift_ui;
        let mut commit = false;
        ui.horizontal(|ui| {
            ui.label(format!("RT shift for {} (min):", self.samples[s].data.name));
            let r = ui.add(dv_rt(&mut v).speed(0.01));
            commit = r.drag_stopped() || r.lost_focus();
            ui.weak("applies to all chromatograms of this file; resets its manual integrations");
        });
        self.samples[s].shift_ui = v;
        if commit && (v - self.samples[s].rt_shift).abs() > 1e-9 {
            self.push_undo("RT shift");
            self.set_shift(s, v, "manual");
        }
    }

    fn ui_results(&mut self, ui: &mut egui::Ui, cur: Option<usize>) {
        let curve = self.mode == Mode::Curve;
        egui::ScrollArea::vertical().id_salt("res").max_height(220.0).show(ui, |ui| {
            egui::Grid::new("restab").striped(true).show(ui, |ui| {
                for h in ["Target", "Start", "End", "Area (a.u.)"] { ui.strong(h); }
                if curve { ui.strong("Conc"); }
                ui.strong("Flags");
                ui.end_row();
                for a in 0..self.analytes.len() {
                    if ui.selectable_label(self.sel_a == a, &self.analytes[a].name).clicked() { self.sel_a = a; }
                    match cur.and_then(|s| self.samples[s].cells.get(a)) {
                        Some(c) => {
                            ui.label(format!("{:.2}", c.start)); ui.label(format!("{:.2}", c.end));
                            ui.label(sci(c.areas.first().copied().unwrap_or(0.0)));
                            if curve { ui.label(cur.and_then(|s| self.conc_of(s, a)).map(sci).unwrap_or_else(|| "-".into())); }
                            ui.label(cur.map(|s| self.flags(s, a)).unwrap_or_default());
                        }
                        None => { for _ in 0..(if curve { 5 } else { 4 }) { ui.label("-"); } }
                    }
                    ui.end_row();
                }
            });
        });
    }

    fn ui_cal(&mut self, ui: &mut egui::Ui) {
        let a = self.sel_a;
        let fit = self.fits.get(a).copied().flatten();
        match fit {
            Some(f) => ui.label(format!("{}: area = {}·conc + {}   R² = {}   n = {}   weighting {}{}", self.analytes[a].name, sci(f.m), sci(f.b), sci(f.r2),
                sci(f.n as f64), WEIGHTS[self.weighting as usize], if f.n < 5 { "   ⚠ fewer than 5 points" } else { "" })),
            None => ui.label("No calibration: mark ≥2 concentration levels as 'Standard' with concentration > 0"),
        };
        let mut toggled = None;
        ui.columns(2, |cols| {
            cols[0].push_id("calt", |ui| {
                egui::Grid::new("caltab").striped(true).show(ui, |ui| {
                    for h in ["Standard", "Conc", "Area", "Accuracy %", "Use"] { ui.strong(h); }
                    ui.end_row();
                    for (si, sm) in self.samples.iter().enumerate() {
                        if sm.kind != Kind::Standard || sm.cells.len() <= a { continue; }
                        let c = &sm.cells[a];
                        ui.label(&sm.data.name); ui.label(sci(sm.conc)); ui.label(sci(c.areas[0]));
                        ui.label(match fit { Some(f) if f.m != 0.0 && sm.conc > 0.0 => sci(100.0 * ((c.areas[0] - f.b) / f.m) / sm.conc), _ => "-".into() });
                        let mut inc = c.include;
                        if ui.checkbox(&mut inc, "").changed() { toggled = Some(si); }
                        ui.end_row();
                    }
                });
            });
            let (mut inc_pts, mut exc_pts) = (vec![], vec![]);
            for sm in self.samples.iter().filter(|s| s.kind == Kind::Standard && s.conc > 0.0 && s.cells.len() > a && s.cells[a].found) {
                let p = [sm.conc, sm.cells[a].areas[0]];
                if sm.cells[a].include { inc_pts.push(p) } else { exc_pts.push(p) }
            }
            sciplot!(Plot::new("cal").height(220.0).x_axis_label(format!("Concentration ({})", self.units)).y_axis_label("Area"))
                .show(&mut cols[1], |pu| {
                    pu.points(Points::new(inc_pts).radius(4.0_f32).name("used"));
                    pu.points(Points::new(exc_pts).radius(4.0_f32).color(Color32::GRAY).name("excluded"));
                    if let Some(f) = fit { pu.line(Line::new(PlotPoints::from(vec![[0.0, f.b], [f.hi * 1.05, f.m * f.hi * 1.05 + f.b]])).name("fit")); }
                });
        });
        if let Some(si) = toggled {
            self.push_undo("calibration point");
            let inc = !self.samples[si].cells[a].include;
            self.samples[si].cells[a].include = inc;
            let msg = format!("file={} analyte={} include={} reason={}", self.samples[si].data.name, self.analytes[a].name, inc,
                if self.reason.trim().is_empty() { "NONE GIVEN" } else { self.reason.trim() });
            self.audit.log("CAL_POINT_TOGGLE", &msg);
            self.refit(a);
        }
    }

    /// Areas or concentrations of every target in every imported file.
    fn ui_matrix(&mut self, ui: &mut egui::Ui) {
        if self.analytes.is_empty() { ui.label("Load a target table to see the results matrix."); return; }
        let curve = self.mode == Mode::Curve;
        let use_conc = curve && self.matrix_conc;
        ui.horizontal(|ui| {
            ui.label("Values:");
            ui.selectable_value(&mut self.matrix_conc, false, "Areas (a.u.)");
            if curve { ui.selectable_value(&mut self.matrix_conc, true, format!("Concentrations ({})", self.units)); }
            if ui.button("Copy as TSV").clicked() {
                let mut t = String::from("file\ttype");
                for a in &self.analytes { t += &format!("\t{}", a.name); }
                t.push('\n');
                for s in 0..self.samples.len() {
                    t += &format!("{}\t{}", self.samples[s].data.name, self.samples[s].kind.label());
                    for a in 0..self.analytes.len() {
                        let c = self.samples[s].cells.get(a);
                        let v = match c { Some(c) if c.found && c.areas[0] > 0.0 => if use_conc { self.conc_of(s, a) } else { Some(c.areas[0]) }, _ => None };
                        t += &match v { Some(v) => format!("\t{:.6e}", v), None => "\tND".to_string() };
                    }
                    t.push('\n');
                }
                ui.ctx().copy_text(t);
                self.status = "Results matrix copied to the clipboard (tab-separated)".into();
            }
            ui.weak("yellow = flagged (see Flags in the chromatogram tab), ND = not detected");
        });
        egui::ScrollArea::both().id_salt("matrix").show(ui, |ui| {
            egui::Grid::new("matrix").striped(true).show(ui, |ui| {
                ui.strong("File"); ui.strong("Type");
                for a in &self.analytes { ui.strong(&a.name); }
                ui.end_row();
                for s in 0..self.samples.len() {
                    ui.label(&self.samples[s].data.name);
                    ui.label(self.samples[s].kind.label());
                    for a in 0..self.analytes.len() {
                        match self.samples[s].cells.get(a) {
                            Some(c) if c.found && c.areas[0] > 0.0 => {
                                let v = if use_conc { self.conc_of(s, a).map(sci).unwrap_or_else(|| "-".into()) } else { sci(c.areas[0]) };
                                let fl = self.flags(s, a);
                                if fl.is_empty() || fl == "manual" { ui.label(v); } else { ui.colored_label(Color32::YELLOW, v).on_hover_text(fl); }
                            }
                            _ => { ui.label("ND"); }
                        }
                    }
                    ui.end_row();
                }
            });
        });
    }

    fn ui_chrom(&mut self, ui: &mut egui::Ui) {
        let cur = self.cur();
        ui.horizontal(|ui| {
            ui.label("View:");
            ui.selectable_value(&mut self.view, View::Xic, "Target XICs");
            ui.selectable_value(&mut self.view, View::Tic, "TIC");
            ui.selectable_value(&mut self.view, View::Bpc, "BPC");
            ui.separator();
            if self.view == View::Xic {
                egui::ComboBox::new("scope", "").width(330.0).selected_text(match self.scope {
                    Scope::Ions => "Selected target, all ions, current file (edit integration)",
                    Scope::AllTargets => "All targets (quant ions), current file",
                    Scope::Files => "Selected target (quant ion), all selected files",
                }).show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.scope, Scope::Ions, "Selected target, all ions, current file (edit integration)");
                    ui.selectable_value(&mut self.scope, Scope::AllTargets, "All targets (quant ions), current file");
                    ui.selectable_value(&mut self.scope, Scope::Files, "Selected target (quant ion), all selected files");
                });
            } else {
                ui.label("Scans:");
                ui.selectable_value(&mut self.tic_level, 1u8, "MS1");
                ui.selectable_value(&mut self.tic_level, 2u8, "MS2");
            }
            ui.checkbox(&mut self.normalize, "Normalize");
            if !self.analytes.is_empty() && self.view == View::Xic && self.scope != Scope::AllTargets {
                ui.separator();
                ui.label("Target:");
                egui::ComboBox::new("an", "").selected_text(self.analytes[self.sel_a].name.clone()).show_ui(ui, |ui| {
                    for i in 0..self.analytes.len() { let n = self.analytes[i].name.clone(); ui.selectable_value(&mut self.sel_a, i, n); }
                });
            }
            ui.label(format!("Current file: {}", cur.map(|s| self.samples[s].data.name.clone()).unwrap_or_else(|| "none".into())));
        });
        let needs_cur = self.view == View::Xic && self.scope != Scope::Files;
        if self.view == View::Xic && self.analytes.is_empty() {
            ui.label("Load a target table (.csv) to see target XICs. TIC and BPC views work without one.");
        } else if needs_cur && cur.is_none() {
            ui.label("No file selected.");
        } else if self.view == View::Xic && self.scope == Scope::Ions {
            if let Some(s) = cur { self.ui_plot(ui, s); }
        } else {
            self.ui_overlay(ui, cur);
        }
        if let Some(s) = cur { self.ui_shift(ui, s); }
        if !self.analytes.is_empty() {
            ui.separator();
            self.ui_results(ui, cur);
            if self.mode == Mode::Curve {
                ui.separator();
                self.ui_cal(ui);
            }
        }
    }

    fn fit_cols(&mut self, avail: f32) {
        let w = (avail - 28.0 - 4.0 * 6.0 - 130.0).max(260.0);
        self.col_w = [w * 0.22, w * 0.10, w * 0.14, w * 0.54];
    }

    fn ui_targets_window(&mut self, ctx: &egui::Context) {
        if !self.show_targets { return; }
        let mut open = true;
        let (mut commit, mut delete): (Option<usize>, Option<usize>) = (None, None);
        let (mut add, mut save, mut save_as, mut open_t) = (false, false, false, false);
        egui::Window::new("Target table (live editable)").open(&mut open).default_size([1000.0, 440.0]).resizable(true).show(ctx, |ui| {
            if self.rows_dirty || self.rows.len() != self.analytes.len() { self.sync_rows(); }
            let avail = ui.available_width();
            if self.col_w[0] < 1.0 { self.fit_cols(avail); }
            ui.horizontal(|ui| {
                if ui.button("Add target").clicked() { add = true; }
                if ui.button("Open…").clicked() { open_t = true; }
                if ui.button("Save").clicked() { save = true; }
                if ui.button("Save as…").clicked() { save_as = true; }
                if ui.button("Fit columns to window").clicked() { self.fit_cols(avail); }
                ui.weak(match &self.table_path { Some(p) => p.display().to_string(), None => "(not saved yet)".into() });
            });
            ui.weak("Drag the column borders to resize; the window scrolls sideways if needed. Changes apply when you leave a field and re-integrate immediately. Precursor m/z is only needed for MS1 quantification.");
            egui::ScrollArea::both().auto_shrink([false, false]).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    ui.add_sized([28.0, 18.0], egui::Label::new(RichText::new("#").strong()));
                    let heads = ["Name", "RT (min)", "Precursor m/z (MS1)", "Product ion m/z (MS2), separated by ;"];
                    for k in 0..4 {
                        ui.add_sized([self.col_w[k], 18.0], egui::Label::new(RichText::new(heads[k]).strong()));
                        let (rect, resp) = ui.allocate_exact_size(egui::vec2(6.0, 18.0), egui::Sense::drag());
                        ui.painter().vline(rect.center().x, rect.y_range(), ui.visuals().widgets.noninteractive.bg_stroke);
                        if resp.hovered() || resp.dragged() { ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeHorizontal); }
                        if resp.dragged() { self.col_w[k] = (self.col_w[k] + resp.drag_delta().x).max(40.0); }
                    }
                });
                for i in 0..self.rows.len() {
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = 0.0;
                        ui.add_sized([28.0, 20.0], egui::Label::new((i + 1).to_string()));
                        for k in 0..4 {
                            let w = self.col_w[k];
                            let r = ui.add_sized([w, 20.0], egui::TextEdit::singleline(&mut self.rows[i][k]));
                            if r.lost_focus() { commit = Some(i); }
                            ui.add_space(6.0);
                        }
                        if parse_row(&self.rows[i]).is_none() { ui.colored_label(Color32::RED, "invalid"); }
                        ui.add_space(8.0);
                        if ui.small_button("Delete").clicked() { delete = Some(i); }
                    });
                }
            });
        });
        self.show_targets = open;
        if let Some(i) = commit { self.apply_row(i); }
        if let Some(i) = delete { self.delete_target(i); }
        if add { self.add_target(); }
        if open_t { self.open_targets_dialog(); }
        if save_as { self.save_targets_as(); }
        if save {
            match self.table_path.clone() { Some(p) => self.write_targets(&p), None => self.save_targets_as() }
        }
    }

    fn menu(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        egui::menu::bar(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.add(egui::Button::new("New project").shortcut_text("Ctrl+N")).clicked() { self.new_project(); ui.close_menu(); }
                if ui.add(egui::Button::new("Open project…").shortcut_text("Ctrl+O")).clicked() { self.open_project_dialog(ctx); ui.close_menu(); }
                if ui.add(egui::Button::new("Save project").shortcut_text("Ctrl+S")).clicked() { self.save_project(false); ui.close_menu(); }
                if ui.button("Save project as…").clicked() { self.save_project(true); ui.close_menu(); }
                ui.separator();
                if ui.button("Import Thermo .raw files…").clicked() {
                    if let Some(files) = rfd::FileDialog::new().set_title("Import Thermo .raw files").add_filter("Thermo RAW", &["raw"]).pick_files() {
                        for f in files { self.start_load(ctx, &f.to_string_lossy()); }
                    }
                    ui.close_menu();
                }
                if ui.button("Open target table…").clicked() { self.open_targets_dialog(); ui.close_menu(); }
                if ui.button("Save target table as…").clicked() { self.save_targets_as(); ui.close_menu(); }
                ui.separator();
                if ui.button("Export analysis results…").clicked() {
                    self.choose_out_path();
                    self.status = match self.export() { Ok(_) => format!("Exported {}", self.out_path), Err(e) => format!("Export failed: {}", e) };
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Exit").clicked() { ctx.send_viewport_cmd(egui::ViewportCommand::Close); }
            });
            ui.menu_button("Edit", |ui| {
                let (nu, nr) = self.history.lens();
                let ul = self.history.peek_undo().map(|s| format!("Undo {}", s.label)).unwrap_or_else(|| "Undo".into());
                let rl = self.history.peek_redo().map(|s| format!("Redo {}", s.label)).unwrap_or_else(|| "Redo".into());
                if ui.add_enabled(nu > 0, egui::Button::new(ul).shortcut_text("Ctrl+Z")).clicked() { self.undo(); ui.close_menu(); }
                if ui.add_enabled(nr > 0, egui::Button::new(rl).shortcut_text("Ctrl+Y")).clicked() { self.redo(); ui.close_menu(); }
                ui.weak(format!("history: {} back / {} forward (max {})", nu, nr, MAX_UNDO));
                ui.separator();
                if ui.button("Edit target table…").clicked() { self.show_targets = true; ui.close_menu(); }
            });
            ui.menu_button("View", |ui| {
                if ui.add(egui::Button::new("Zoom in").shortcut_text("Ctrl++")).clicked() { self.zoom = (self.zoom + 0.1).min(3.0); }
                if ui.add(egui::Button::new("Zoom out").shortcut_text("Ctrl+-")).clicked() { self.zoom = (self.zoom - 0.1).max(0.5); }
                if ui.add(egui::Button::new("Reset zoom").shortcut_text("Ctrl+0")).clicked() { self.zoom = 1.0; }
                ui.separator();
                if ui.selectable_label(self.tab == Tab::Chrom, "Chromatograms").clicked() { self.tab = Tab::Chrom; ui.close_menu(); }
                if ui.selectable_label(self.tab == Tab::Matrix, "Results matrix (all files)").clicked() { self.tab = Tab::Matrix; ui.close_menu(); }
                ui.separator();
                ui.checkbox(&mut self.show_spectrum, "MS2 spectrum panel");
                if ui.button("Audit trail").clicked() { self.show_audit = true; ui.close_menu(); }
            });
            ui.menu_button("Analysis", |ui| {
                if ui.button("Re-integrate all (apply settings)").clicked() { self.apply_settings(); ui.close_menu(); }
                if ui.button("Align retention times now").clicked() { self.align_all(); ui.close_menu(); }
                if ui.button("Clear retention-time alignment").clicked() { self.clear_alignment(); ui.close_menu(); }
                ui.separator();
                if ui.selectable_label(self.mode == Mode::Curve, "Quantify with standard curve").clicked() { self.mode = Mode::Curve; ui.close_menu(); }
                if ui.selectable_label(self.mode == Mode::Areas, "Report integrated areas only").clicked() { self.mode = Mode::Areas; ui.close_menu(); }
                ui.separator();
                if ui.selectable_label(self.basis == Basis::Ms2, "Quantify on MS2 product-ion area").clicked() { self.set_basis(Basis::Ms2); ui.close_menu(); }
                if ui.selectable_label(self.basis == Basis::Ms1, "Quantify on MS1 precursor area").clicked() { self.set_basis(Basis::Ms1); ui.close_menu(); }
                ui.separator();
                if ui.button("Verify audit trail integrity").clicked() {
                    self.status = match audit::verify(&self.audit.path) {
                        Ok(n) => format!("Audit trail intact: {} entries, hash chain verified", n),
                        Err(e) => format!("AUDIT TRAIL PROBLEM: {}", e),
                    };
                    ui.close_menu();
                }
            });
            ui.menu_button("Help", |ui| {
                if ui.button("README (offline)").clicked() { self.show_readme = true; ui.close_menu(); }
                if ui.button("Credits & licenses").clicked() { self.show_credits = true; ui.close_menu(); }
                ui.separator();
                ui.label(format!("LC-MS/MS Quant v{}", VERSION));
            });
        });
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        ctx.options_mut(|o| o.zoom_with_keyboard = false);
        let (zin, zout, zreset) = ctx.input_mut(|i| (
            i.consume_key(Modifiers::COMMAND, Key::Plus) || i.consume_key(Modifiers::COMMAND, Key::Equals),
            i.consume_key(Modifiers::COMMAND, Key::Minus),
            i.consume_key(Modifiers::COMMAND, Key::Num0),
        ));
        if zin { self.zoom = (self.zoom + 0.1).min(3.0); }
        if zout { self.zoom = (self.zoom - 0.1).max(0.5); }
        if zreset { self.zoom = 1.0; }
        if (self.zoom - ctx.zoom_factor()).abs() > 1e-4 { ctx.set_zoom_factor(self.zoom); self.write_settings(); }
        if ctx.wants_keyboard_input() { return; }
        let cs = Modifiers::COMMAND.plus(Modifiers::SHIFT);
        let (redo, undo, save, open, new) = ctx.input_mut(|i| (
            i.consume_key(cs, Key::Z) || i.consume_key(Modifiers::COMMAND, Key::Y),
            i.consume_key(Modifiers::COMMAND, Key::Z),
            i.consume_key(Modifiers::COMMAND, Key::S),
            i.consume_key(Modifiers::COMMAND, Key::O),
            i.consume_key(Modifiers::COMMAND, Key::N),
        ));
        if redo { self.redo(); }
        if undo { self.undo(); }
        if save { self.save_project(false); }
        if open { self.open_project_dialog(ctx); }
        if new { self.new_project(); }
    }
}

impl App {
    /// m/z values traced for target `a`: the product ions (MS2 basis) or the precursor (MS1 basis).
    fn trace_mz(&self, a: usize) -> Vec<f64> {
        match self.basis {
            Basis::Ms2 => self.analytes[a].ions.clone(),
            Basis::Ms1 => vec![self.analytes[a].precursor.unwrap_or(0.0)],
        }
    }

    fn set_basis(&mut self, b: Basis) {
        if self.basis == b { return; }
        self.push_undo("quantification signal");
        self.basis = b;
        self.audit.log("QUANT_BASIS", if b == Basis::Ms1 { "MS1 precursor area" } else { "MS2 product-ion area" });
        self.recompute_all();
        if b == Basis::Ms1 && self.analytes.iter().any(|a| a.precursor.is_none()) {
            self.status = "Some targets have no precursor m/z (add p:<m/z> in the table); they cannot be quantified on the MS1 signal.".into();
        }
    }

    /// MS2 scan (index into `spectra`) closest in time to `rt` (aligned minutes) in file `f`.
    fn nearest_ms2(&self, f: usize, rt: f64) -> Option<usize> {
        let sm = self.samples.get(f)?;
        let sp = &sm.data.spectra;
        if sp.is_empty() { return None; }
        let t = rt - sm.rt_shift;
        let j = sp.partition_point(|x| x.rt < t);
        let (c0, c1) = (j.saturating_sub(1), j.min(sp.len() - 1));
        Some(if (sp[c1].rt - t).abs() < (sp[c0].rt - t).abs() { c1 } else { c0 })
    }

    /// MS2 scan inside the integration window with the strongest quantifier product ion.
    fn apex_scan(&self, s: usize, a: usize) -> Option<usize> {
        let sm = self.samples.get(s)?;
        let c = sm.cells.get(a)?;
        if !c.found { return None; }
        let ion = *self.analytes.get(a)?.ions.first()?;
        let tol = self.tol_da(ion);
        let (t0, t1) = (c.start - sm.rt_shift, c.end - sm.rt_shift);
        let sp = &sm.data.spectra;
        let (mut best, mut top) = (None, 0.0f64);
        for j in sp.partition_point(|x| x.rt < t0)..sp.len() {
            if sp[j].rt > t1 { break; }
            let v = intensity_near(&sp[j].mz, &sp[j].inten, ion, tol);
            if v > top { top = v; best = Some(j); }
        }
        best
    }

    fn pin_scan(&mut self, f: usize, rt: f64) {
        if let Some(i) = self.nearest_ms2(f, rt) { self.pinned_scan = Some((f, i)); }
    }

    /// Which MS2 scan the spectrum panel shows: hovered > pinned > apex of the selected target in the current file.
    fn spectrum_pick(&self) -> Option<(usize, usize, &'static str)> {
        if self.follow_hover {
            if let Some((f, rt)) = self.hover { if let Some(i) = self.nearest_ms2(f, rt) { return Some((f, i, "hovered")); } }
        }
        if let Some((f, i)) = self.pinned_scan {
            if f < self.samples.len() && i < self.samples[f].data.spectra.len() { return Some((f, i, "pinned")); }
        }
        let s = self.cur()?;
        if self.sel_a < self.analytes.len() { return self.apex_scan(s, self.sel_a).map(|i| (s, i, "apex of the selected target")); }
        None
    }

    /// Second graph: MS2 spectrum of the selected target (stems; target product ions highlighted).
    fn ui_spectrum(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.strong("MS2 spectrum");
            ui.checkbox(&mut self.follow_hover, "follow chromatogram hover");
            ui.checkbox(&mut self.spec_autofit, "auto-fit view");
            if self.pinned_scan.is_some() && ui.button("Unpin scan").clicked() { self.pinned_scan = None; }
            ui.weak("click the chromatogram to pin the nearest MS2 scan");
        });
        let (f, i, why) = match self.spectrum_pick() {
            Some(x) => x,
            None => { ui.label("No MS2 scan to show: select a file and a target with a detected peak, or hover / click the chromatogram."); return; }
        };
        let sm = &self.samples[f];
        let sp = &sm.data.spectra[i];
        let norm = self.normalize;
        let mx = sp.inten.iter().cloned().fold(0.0f32, f32::max) as f64;
        let scale = if norm && mx > 0.0 { 1.0 / mx } else { 1.0 };
        ui.label(format!("{} · MS2 scan {} · RT {:.2} min · {} peaks · shown: {}", sm.data.name, sp.scan, sp.rt + sm.rt_shift, sp.mz.len(), why));
        let mut all: Vec<[f64; 2]> = Vec::with_capacity(sp.mz.len() * 3);
        for (m, y) in sp.mz.iter().zip(sp.inten.iter()) { let y = *y as f64 * scale; all.push([*m, 0.0]); all.push([*m, y]); all.push([*m, 0.0]); }
        let mut hl: Vec<(usize, f64, Vec<[f64; 2]>)> = vec![];
        if self.sel_a < self.analytes.len() {
            for (k, &ion) in self.analytes[self.sel_a].ions.iter().enumerate() {
                let tol = self.tol_da(ion);
                let mut v = vec![];
                for (m, y) in sp.mz.iter().zip(sp.inten.iter()) {
                    if (*m - ion).abs() <= tol { let y = *y as f64 * scale; v.push([*m, 0.0]); v.push([*m, y]); v.push([*m, 0.0]); }
                }
                hl.push((k, ion, v));
            }
        }
        let mut plot = sciplot!(Plot::new("spec").height(ui.available_height()).legend(Legend::default()).x_axis_label("m/z")
            .y_axis_label(if norm { "Relative intensity" } else { "Intensity" }));
        if self.spec_autofit { plot = plot.reset(); }
        plot.show(ui, |pu| {
            pu.line(Line::new(PlotPoints::from(all)).color(Color32::GRAY).name("peaks"));
            for (k, ion, v) in hl {
                pu.vline(VLine::new(ion).color(pal(k).gamma_multiply(0.5)));
                if !v.is_empty() { pu.line(Line::new(PlotPoints::from(v)).color(pal(k)).width(2.5_f32).name(format!("target ion m/z {}", sci_prec(ion, 9)))); }
            }
        });
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _f: &mut eframe::Frame) {
        while let Ok(m) = self.rx.try_recv() {
            match m {
                Msg::Loaded(d, info) => self.add_sample(d, info),
                Msg::Err(e) => self.status = format!("Load failed: {}", e),
                Msg::Finished => self.loading = self.loading.saturating_sub(1),
            }
        }
        if self.loading == 0 { if let Some(p) = self.pending.take() { self.finish_project(p); } }
        for f in ctx.input(|i| i.raw.dropped_files.clone()) {
            if let Some(p) = f.path { self.load_path(ctx, &p.to_string_lossy()); }
        }
        self.shortcuts(ctx);
        self.fix_selection();
        if self.sel_a >= self.analytes.len() { self.sel_a = 0; }
        if self.align_ref >= self.samples.len() { self.align_ref = 0; }

        let title = format!("LC-MS/MS Quant - {}", self.project_path.as_ref().and_then(|p| p.file_name()).map(|n| n.to_string_lossy().to_string()).unwrap_or_else(|| "untitled project".into()));
        if title != self.title { ctx.send_viewport_cmd(egui::ViewportCommand::Title(title.clone())); self.title = title; }

        egui::TopBottomPanel::top("menu").show(ctx, |ui| { self.menu(ui, ctx); });

        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("User:");
                let r = ui.add(egui::TextEdit::singleline(&mut self.audit.user).desired_width(110.0));
                if r.lost_focus() {
                    if self.audit.user.trim().is_empty() { self.audit.user = "N/A".into(); }
                    if self.audit.user != self.user_logged {
                        let msg = format!("{} -> {}", self.user_logged, self.audit.user);
                        self.audit.log("USER_CHANGED", &msg);
                        self.user_logged = self.audit.user.clone();
                        self.write_settings();
                    }
                }
                ui.separator();
                ui.label("File path:");
                ui.add(egui::TextEdit::singleline(&mut self.path_in).desired_width(260.0));
                if ui.button("Load").clicked() { let p = self.path_in.trim().trim_matches('"').to_string(); self.load_path(ctx, &p); }
                ui.separator();
                ui.label("Results file:");
                let r = ui.add(egui::TextEdit::singleline(&mut self.out_path).desired_width(220.0)).on_hover_text("Click to choose the output file in the file browser");
                if r.clicked() { self.choose_out_path(); }
                if ui.button("Export CSV").clicked() {
                    self.status = match self.export() { Ok(_) => format!("Exported {}", self.out_path), Err(e) => format!("Export failed: {}", e) };
                }
                ui.separator();
                ui.label("ThermoRawFileParser.exe:");
                let r = ui.add(egui::TextEdit::singleline(&mut self.trfp_path).desired_width(260.0)
                    .hint_text(r"C:\Tools\ThermoRawFileParser\ThermoRawFileParser.exe"));
                if r.lost_focus() { self.save_trfp(); }
                if !self.trfp_path.trim().is_empty() {
                    if Path::new(self.trfp_path.trim().trim_matches('"')).is_file() { ui.colored_label(Color32::GREEN, "found"); }
                    else { ui.colored_label(Color32::RED, "not found"); }
                }
                if self.loading > 0 { ui.spinner(); }
            });
            ui.horizontal(|ui| {
                let mut ch = false;
                egui::ComboBox::new("tolu", "").selected_text(if self.tol_ppm { "ppm" } else { "Da" }).show_ui(ui, |ui| {
                    ch |= ui.selectable_value(&mut self.tol_ppm, true, "ppm").changed();
                    ch |= ui.selectable_value(&mut self.tol_ppm, false, "Da").changed();
                });
                ui.label("m/z tol ±");
                ch |= ui.add(dv(&mut self.tol).speed(0.1).range(0.001..=1000.0)).changed();
                ui.label("RT window ± (min)");
                ch |= ui.add(dv_rt(&mut self.rt_win).speed(0.01).range(0.01..=60.0)).changed();
                ui.label("Weighting");
                egui::ComboBox::new("wt", "").selected_text(WEIGHTS[self.weighting as usize]).show_ui(ui, |ui| {
                    for (i, w) in WEIGHTS.iter().enumerate() { ch |= ui.selectable_value(&mut self.weighting, i as u8, *w).changed(); }
                });
                ch |= ui.checkbox(&mut self.baseline, "Linear baseline subtraction").changed();
                ui.label("Ion ratio tol ±%");
                ch |= ui.add(dv(&mut self.ratio_tol).range(1.0..=100.0)).changed();
                ui.label("Units");
                ui.add(egui::TextEdit::singleline(&mut self.units).desired_width(60.0));
                if ch { self.settings_dirty = true; }
                if ui.add_enabled(self.settings_dirty, egui::Button::new("Apply & re-integrate")).clicked() { self.apply_settings(); }
                if self.settings_dirty { ui.colored_label(Color32::YELLOW, "settings changed - not applied"); }
            });
            ui.horizontal(|ui| {
                ui.label("Quantification:");
                let before = self.mode;
                egui::ComboBox::new("mode", "").selected_text(if self.mode == Mode::Areas { "Integrated areas only (a.u.)" } else { "Standard curve" }).show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.mode, Mode::Areas, "Integrated areas only (a.u.)");
                    ui.selectable_value(&mut self.mode, Mode::Curve, "Standard curve");
                });
                if self.mode != before { self.audit.log("QUANT_MODE", if self.mode == Mode::Areas { "integrated areas only" } else { "standard curve" }); }
                ui.label("Signal:");
                let mut nb = self.basis;
                egui::ComboBox::new("basis", "").selected_text(if self.basis == Basis::Ms2 { "MS2 product-ion area" } else { "MS1 precursor area" }).show_ui(ui, |ui| {
                    ui.selectable_value(&mut nb, Basis::Ms2, "MS2 product-ion area");
                    ui.selectable_value(&mut nb, Basis::Ms1, "MS1 precursor area");
                });
                if nb != self.basis { self.set_basis(nb); }
                ui.separator();
                if ui.checkbox(&mut self.align_on, "Align RT on import").changed() { self.audit.log("ALIGN_SETTING", &format!("align_on_import={}", self.align_on)); }
                ui.label("max shift ±(min)");
                ui.add(dv_rt(&mut self.align_max).speed(0.05).range(0.1..=10.0));
                ui.label("reference:");
                egui::ComboBox::new("aref", "")
                    .selected_text(self.samples.get(self.align_ref).map(|s| s.data.name.clone()).unwrap_or_else(|| "-".into()))
                    .show_ui(ui, |ui| { for i in 0..self.samples.len() { ui.selectable_value(&mut self.align_ref, i, &self.samples[i].data.name); } });
                if ui.button("Align all now").clicked() { self.align_all(); }
                if ui.button("Clear alignment").clicked() { self.clear_alignment(); }
                ui.separator();
                if ui.button("Targets…").clicked() { self.show_targets = true; }
                ui.weak(format!("zoom {}%", (self.zoom * 100.0).round()));
            });
        });
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| { ui.label(&self.status); });

        egui::SidePanel::left("samples").min_width(400.0).show(ctx, |ui| {
            ui.heading("Samples");
            ui.horizontal(|ui| {
                if ui.small_button("Select all").clicked() { for s in self.samples.iter_mut() { s.show = true; } }
                if ui.small_button("Select none").clicked() { for s in self.samples.iter_mut() { s.show = false; } }
            });
            ui.weak("Click = select · Ctrl+click = add/remove · Shift+click = range · right-click = menu");
            let before: Vec<(Kind, f64)> = self.samples.iter().map(|s| (s.kind, s.conc)).collect();
            let mut changed = false;
            let mut click: Option<(usize, bool, bool)> = None;
            let mut remove: Option<usize> = None;
            let mut set_ref: Option<usize> = None;
            egui::ScrollArea::vertical().show(ui, |ui| {
                for (i, s) in self.samples.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut s.show, "");
                        let label = if i == self.sel_s && s.show { RichText::new(&s.data.name).strong() } else { RichText::new(&s.data.name) };
                        let resp = ui.selectable_label(s.show, label);
                        if resp.clicked() { let m = ui.input(|inp| inp.modifiers); click = Some((i, m.command, m.shift)); }
                        resp.context_menu(|ui| {
                            if ui.button("Remove from project (file stays on disk)").clicked() { remove = Some(i); ui.close_menu(); }
                            if ui.button("Select only this file").clicked() { click = Some((i, false, false)); ui.close_menu(); }
                            if ui.button("Use as alignment reference").clicked() { set_ref = Some(i); ui.close_menu(); }
                        });
                        egui::ComboBox::new(("kind", i), "").selected_text(s.kind.label()).show_ui(ui, |ui| {
                            for k in Kind::ALL { changed |= ui.selectable_value(&mut s.kind, k, k.label()).changed(); }
                        });
                        if self.mode == Mode::Curve && s.kind == Kind::Standard { changed |= ui.add(dv(&mut s.conc).speed(0.1).range(0.0..=1e9)).changed(); }
                        if s.rt_shift != 0.0 { ui.weak(format!("Δ{:+.2} min", s.rt_shift)); }
                    });
                }
            });
            if changed {
                let dragging = ctx.input(|i| i.pointer.any_down());
                if !(dragging && self.last_undo_label == "sample assignment") {
                    let mut snap = self.snapshot("sample assignment");
                    for (st, (k, c)) in snap.samples.iter_mut().zip(before.iter()) { st.kind = *k; st.conc = *c; }
                    self.history.record(snap);
                    self.last_undo_label = "sample assignment".into();
                }
                for a in 0..self.analytes.len() { self.refit(a); }
                self.audit.log("SAMPLE_ASSIGNMENT", &self.samples.iter().map(|s| format!("{}={}@{}", s.data.name, s.kind.label(), s.conc)).collect::<Vec<_>>().join(" ; "));
            }
            if let Some((i, ctrl, shift)) = click { self.apply_click(i, ctrl, shift); }
            if let Some(i) = set_ref { self.align_ref = i; self.audit.log("ALIGN_REFERENCE", &self.samples[i].data.name); }
            if let Some(i) = remove { self.remove_sample(i); }
        });

        if self.show_spectrum && self.tab == Tab::Chrom && !self.samples.is_empty() {
            egui::TopBottomPanel::bottom("spectrum").resizable(true).default_height(270.0).min_height(120.0).show(ctx, |ui| { self.ui_spectrum(ui); });
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            if self.samples.is_empty() {
                ui.heading("Drop Thermo .raw files (and optionally a .csv target table: name, RT|n/a, ion m/z values) onto this window.");
                ui.label("Or use File → Import / Open. Files are read via ThermoRawFileParser and kept in memory only; no converted files are written to disk.");
                return;
            }
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.tab, Tab::Chrom, "Chromatograms");
                ui.selectable_value(&mut self.tab, Tab::Matrix, "Results matrix (all files)");
            });
            ui.separator();
            let prev_hover = self.hover;
            self.hover = None;
            match self.tab { Tab::Chrom => self.ui_chrom(ui), Tab::Matrix => self.ui_matrix(ui) }
            if self.hover != prev_hover { ui.ctx().request_repaint(); }
        });

        self.ui_targets_window(ctx);
        if self.show_credits {
            let mut open = true;
            egui::Window::new("Credits & third-party software").open(&mut open).default_width(640.0).show(ctx, |ui| { ui.label(CREDITS); });
            self.show_credits = open;
        }
        if self.show_readme {
            let mut open = true;
            egui::Window::new("README (offline)").open(&mut open).default_size([760.0, 560.0]).show(ctx, |ui| {
                egui::ScrollArea::both().show(ui, |ui| { ui.monospace(README); });
            });
            self.show_readme = open;
        }
        if self.show_audit {
            let mut open = true;
            egui::Window::new("Audit trail (hash-chained)").open(&mut open).default_size([900.0, 400.0]).show(ctx, |ui| {
                egui::ScrollArea::both().stick_to_bottom(true).show(ui, |ui| {
                    for l in self.audit.recent.iter().rev().take(300).rev() { ui.monospace(l); }
                });
            });
            self.show_audit = open;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testkit;

    fn app() -> App {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static N: AtomicUsize = AtomicUsize::new(0);
        let p = std::env::temp_dir().join(format!("lcmsq_app_test_{}_{}.log", std::process::id(), N.fetch_add(1, Ordering::SeqCst)));
        let _ = std::fs::remove_file(&p);
        let mut a = App::new_with_audit(&p.to_string_lossy());
        a.trfp_path.clear();
        a
    }
    fn info(n: usize) -> LoadInfo { LoadInfo { path: format!("synthetic{}.raw", n), raw_sha: "0".repeat(64), detail: String::new() } }
    fn add(a: &mut App, n: usize, offset: f64, scale: f64) {
        let d = mzml::parse_bytes(&format!("run{}.raw", n), &testkit::synthetic_mzml(offset, scale)).expect("synthetic mzML parses");
        a.add_sample(d, info(n));
    }
    fn targets(a: &mut App) {
        a.analytes = parse_table("name,rt,ions\nSynth,5.0,150;250\n").unwrap();
        a.recompute_all();
    }
    fn shown(a: &App) -> Vec<usize> { a.shown() }

    #[test]
    fn explorer_style_selection() {
        let mut a = app();
        for i in 0..4 { add(&mut a, i, 0.0, 1.0); }
        a.apply_click(1, false, false);
        assert_eq!(shown(&a), vec![1]);
        a.apply_click(3, true, false);
        assert_eq!(shown(&a), vec![1, 3]);
        a.apply_click(0, false, true); // range from anchor(3) to 0
        assert_eq!(shown(&a), vec![0, 1, 2, 3]);
        a.apply_click(2, true, false); // toggle off
        assert_eq!(shown(&a), vec![0, 1, 3]);
        assert!(a.cur().is_some());
    }

    #[test]
    fn no_selection_means_no_current_file() {
        let mut a = app();
        for i in 0..2 { add(&mut a, i, 0.0, 1.0); }
        for s in a.samples.iter_mut() { s.show = false; }
        a.fix_selection();
        assert!(a.cur().is_none());
        assert!(a.shown().is_empty());
    }

    #[test]
    fn remove_file_keeps_disk_untouched_and_is_undoable() {
        let mut a = app();
        for i in 0..3 { add(&mut a, i, 0.0, 1.0); }
        a.remove_sample(1);
        assert_eq!(a.samples.len(), 2);
        assert_eq!(a.samples[1].data.name, "run2.raw");
        a.undo();
        assert_eq!(a.samples.len(), 3);
        assert_eq!(a.samples[1].data.name, "run1.raw");
        a.redo();
        assert_eq!(a.samples.len(), 2);
    }

    #[test]
    fn undo_history_is_limited_to_ten_steps() {
        let mut a = app();
        targets(&mut a);
        for _ in 0..12 { a.add_target(); }
        assert_eq!(a.analytes.len(), 13);
        for _ in 0..20 { a.undo(); }
        assert_eq!(a.analytes.len(), 3, "only 10 of the 12 additions can be undone");
        for _ in 0..20 { a.redo(); }
        assert_eq!(a.analytes.len(), 13);
    }

    #[test]
    fn live_target_editing_recomputes_areas() {
        let mut a = app();
        add(&mut a, 0, 0.0, 1.0);
        targets(&mut a);
        let first = a.samples[0].cells[0].areas[0];
        assert!(first > 0.0);
        a.sync_rows();
        a.rows[0][3] = "400".into(); // constant 10-count ion: different area
        a.apply_row(0);
        assert_eq!(a.analytes[0].ions, vec![400.0]);
        assert!((a.samples[0].cells[0].areas[0] - first).abs() > 1.0);
        a.rows[0][3] = "".into(); // invalid: must not be applied
        a.apply_row(0);
        assert_eq!(a.analytes[0].ions, vec![400.0]);
        a.add_target();
        assert_eq!(a.analytes.len(), 2);
        assert_eq!(a.samples[0].cells.len(), 2);
        a.delete_target(0);
        assert_eq!(a.analytes.len(), 1);
        assert_eq!(a.samples[0].cells.len(), 1);
        let csv_text = a.targets_csv();
        assert!(parse_table(&csv_text).is_ok());
    }

    #[test]
    fn standard_curve_quantifies_a_known_sample() {
        let mut a = app();
        targets(&mut a);
        for (i, s) in [1.0, 2.0, 4.0, 8.0].iter().enumerate() {
            add(&mut a, i, 0.0, *s);
            a.samples[i].kind = Kind::Standard;
            a.samples[i].conc = *s;
        }
        add(&mut a, 9, 0.0, 3.0);
        for k in 0..a.analytes.len() { a.refit(k); }
        let f = a.fits[0].expect("calibration fitted");
        assert!(f.r2 > 0.9999, "r2 = {}", f.r2);
        let c = a.conc_of(4, 0).expect("unknown quantified");
        assert!((c - 3.0).abs() / 3.0 < 0.02, "got {}", c);
        a.mode = Mode::Areas;
        assert!(a.conc_of(4, 0).is_none(), "areas-only mode reports no concentrations");
        assert!(a.samples[4].cells[0].areas[0] > 0.0);
    }

    #[test]
    fn ion_ratio_and_area_are_close_to_analytical_values() {
        let mut a = app();
        targets(&mut a);
        add(&mut a, 0, 0.0, 2.0);
        assert!(a.samples[0].cells[0].found);
        let expect = testkit::expected_area_seconds(2.0);
        let with_baseline = a.samples[0].cells[0].areas[0];
        // Without baseline subtraction only the tails beyond the automatic window (5 % of the apex) are missing.
        a.baseline = false;
        a.recompute_all();
        let c = &a.samples[0].cells[0];
        assert!((c.areas[0] - expect).abs() / expect < 0.03, "area {} vs {}", c.areas[0], expect);
        assert!((c.areas[1] / c.areas[0] - 0.5).abs() < 0.01);
        // With baseline subtraction the 5 % pedestal at the window edges is removed as well: lower, but not by much.
        assert!(with_baseline < c.areas[0] && with_baseline > 0.8 * expect, "{} vs {}", with_baseline, expect);
    }

    #[test]
    fn alignment_on_import_corrects_a_shifted_run() {
        let mut a = app();
        a.align_on = true;
        targets(&mut a);
        add(&mut a, 0, 0.0, 1.0);
        add(&mut a, 1, 0.2, 1.0);
        assert!((a.samples[1].rt_shift + 0.2).abs() < 0.02, "shift {}", a.samples[1].rt_shift);
        // after alignment both runs integrate to the same area
        let (x, y) = (a.samples[0].cells[0].areas[0], a.samples[1].cells[0].areas[0]);
        assert!((x - y).abs() / x < 0.03);
    }

    #[test]
    fn project_round_trip_keeps_manual_integration() {
        let mut a = app();
        targets(&mut a);
        add(&mut a, 0, 0.0, 1.0);
        a.samples[0].kind = Kind::Standard;
        a.samples[0].conc = 5.0;
        a.apply_edit(0, 0, 4.9, 5.1);
        a.commit(0, 0);
        let pf = a.project_file();
        let back = ProjectFile::parse(&pf.to_text()).unwrap();
        assert_eq!(pf, back);
        assert!(back.samples[0].cells[0].manual);
        let mut b = app();
        b.apply_params(&back.params);
        assert_eq!(b.weighting, a.weighting);
        assert_eq!(b.tol, a.tol);
    }

    #[test]
    fn manual_integration_is_audited_with_old_and_new_values() {
        let mut a = app();
        targets(&mut a);
        add(&mut a, 0, 0.0, 1.0);
        a.reason = "test".into();
        a.apply_edit(0, 0, 4.8, 5.2);
        a.commit(0, 0);
        assert!(a.audit.recent.iter().any(|l| l.contains("MANUAL_INTEGRATION") && l.contains("reason=test")));
        assert!(audit::verify(&a.audit.path).is_ok());
    }

    #[test]
    fn ms1_basis_quantifies_the_precursor_signal() {
        let mut a = app();
        a.analytes = parse_table("name,rt,ions\nSynth,5.0,150;250;p:500\n").unwrap();
        a.recompute_all();
        add(&mut a, 0, 0.0, 1.0);
        let ms2 = a.samples[0].cells[0].areas[0];
        assert_eq!(a.samples[0].xics[0].len(), 2);
        a.set_basis(Basis::Ms1);
        assert_eq!(a.samples[0].xics[0].len(), 1, "one precursor trace");
        let ms1 = a.samples[0].cells[0].areas[0];
        let expect = 50000.0 * 0.1 * (2.0 * std::f64::consts::PI).sqrt() * 60.0;
        assert!((ms1 - expect).abs() / expect < 0.15, "MS1 area {} vs {}", ms1, expect);
        assert!((ms1 - ms2).abs() > 1.0, "the two signals are different measurements");
        a.undo();
        assert!(a.basis == Basis::Ms2);
        assert_eq!(a.samples[0].xics[0].len(), 2);
    }

    #[test]
    fn target_without_precursor_gives_no_peak_on_ms1_basis() {
        let mut a = app();
        targets(&mut a);
        add(&mut a, 0, 0.0, 1.0);
        a.set_basis(Basis::Ms1);
        assert!(!a.samples[0].cells[0].found);
        assert!(a.status.contains("precursor"));
    }

    #[test]
    fn hover_picks_the_nearest_trace_and_point() {
        let mk = |name: &str, y: f64, ci: usize| Series {
            name: name.into(), data: (0..100).map(|i| [i as f64 * 0.1, y]).collect(), ci, cur: false, file: 0, level: 2,
            fname: "f".into(), desc: name.into(), scale: 1.0,
        };
        let s = vec![mk("low", 10.0, 0), mk("high", 1000.0, 1)];
        let h = nearest_hover(&s, 5.04, 900.0, 10.0).unwrap();
        assert_eq!((h.desc.as_str(), h.idx), ("high", 50));
        assert_eq!(nearest_hover(&s, 5.04, 20.0, 10.0).unwrap().desc, "low");
        assert!(nearest_hover(&s, 50.0, 10.0, 10.0).is_none(), "cursor far outside the data");
    }

    #[test]
    fn spectrum_panel_scan_selection() {
        let mut a = app();
        targets(&mut a);
        add(&mut a, 0, 0.0, 1.0);
        let i = a.apex_scan(0, 0).expect("apex scan");
        assert!((a.samples[0].data.spectra[i].rt - 5.0).abs() < 0.02);
        assert!(matches!(a.spectrum_pick(), Some((0, _, "apex of the selected target"))));
        let j = a.nearest_ms2(0, 5.0).unwrap();
        assert!((a.samples[0].data.spectra[j].rt - 5.0).abs() <= 0.005 + 1e-9);
        a.pin_scan(0, 7.0);
        assert!(matches!(a.spectrum_pick(), Some((0, _, "pinned"))));
        a.hover = Some((0, 2.0));
        assert!(matches!(a.spectrum_pick(), Some((0, _, "hovered"))));
        a.follow_hover = false;
        assert!(matches!(a.spectrum_pick(), Some((0, _, "pinned"))));
    }
}
