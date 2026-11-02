#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]
use eframe::egui;
use lcms_quant::{app::App, icon};

fn main() -> eframe::Result<()> {
    let mut vp = egui::ViewportBuilder::default().with_inner_size([1500.0, 900.0]).with_title("LC-MS/MS Quant");
    if let Some((rgba, width, height)) = icon::decode_ico(include_bytes!("../img/lcms-quant-ico.ico")) {
        vp = vp.with_icon(egui::IconData { rgba, width, height });
    }
    let opts = eframe::NativeOptions { viewport: vp, ..Default::default() };
    eframe::run_native("LC-MS/MS Quant", opts, Box::new(|_cc| Ok(Box::new(App::new()))))
}
