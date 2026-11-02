# LC-MS/MS Quant

Desktop app (Windows, macOS) for **targeted quantification of LC-MS/MS data from Thermo `.raw` files**.
You provide a target table (name, retention time or `n/a`, diagnostic fragment-ion m/z values); the program extracts the
fragment-ion chromatograms from the MS2 scans, integrates the peaks (boundaries adjustable by hand, values update live) and
reports either **integrated areas (arbitrary units)** or **concentrations from a standard curve** measured in the same data set.

> Intended for method development, research and routine targeted analysis. The software provides data-integrity
> features (hash-chained audit trail, file hashes, reasons for manual changes) but is **not validated or certified** for
> GxP / 21 CFR Part 11 use. Validation of the software and of the converter version is the user's responsibility.

## Features
- Direct import of Thermo `.raw` files (via ThermoRawFileParser); converted data is kept **in memory only**
- Target table: load, edit live inside the program, save under any name
- Manual integration by dragging boundaries or typing values; baseline subtraction; ion-ratio checks
- Quantification on the **MS2 product-ion area** or the **MS1 precursor area**; calibration (linear, weighting none / 1/x / 1/x²), areas-only mode, results matrix for all files
- Retention-time alignment (global shift by TIC cross-correlation), TIC / BPC / XIC inspection of one or many files
- Hover read-out on every chromatogram (m/z, RT, scan number, intensity of the nearest point) and a linked **MS2 spectrum panel** (hovered / pinned / apex scan, target ions highlighted)
- Explorer-style multi-select (Ctrl / Shift), right-click to remove a file from the project (never deletes it on disk)
- Projects (`.quant`), undo/redo (10 steps), scientific-notation display, zoom with Ctrl +/-/0, user name in the audit trail

## Build
Requirements: [Rust](https://rustup.rs) (stable); on Windows additionally the *Visual Studio Build Tools* with the
"Desktop development with C++" workload; on macOS the Xcode command-line tools.

```powershell
git clone <your-repo-url>
cd lcms-quant
cargo run --release          # or: cargo build --release  -> target/release/lcms-quant(.exe)
```

### Thermo reader (needed at run time)
1. Download a release of [ThermoRawFileParser](https://github.com/CompOmics/ThermoRawFileParser/releases) and unzip it to
   `tools/ThermoRawFileParser<anything>/` (or enter the path to `ThermoRawFileParser.exe` in the top bar).
2. The framework-dependent builds need the **.NET 8 runtime**; the self-contained builds do not.
   Test: `dotnet ThermoRawFileParser.dll -v` or `ThermoRawFileParser.exe -v`.

## Use
1. Drop a target table `.csv` and `.raw` files onto the window (or use the File menu). Table format: `name,rt_min,ions` – `n/a` for an unknown RT,
   MS2 product ions separated by `;` (the first is the quantifier) and, optionally, the MS1 precursor as `p:286.1438`, e.g.
   `Morphine,n/a,286.14;201.07;165.05;p:286.1438`. The precursor is only needed for MS1-based quantification. Mark calibration files as **Standard** and enter their concentrations.
2. Choose the signal (MS2 product ions or MS1 precursor) in the top bar. Check/adjust the integrated peaks; drag the green/red lines.
   Give a reason for manual changes. Hover the chromatogram to see m/z, RT and scan; click it to pin the nearest MS2 scan in the spectrum panel.
3. Export the results (`Export CSV`), save the project (`Ctrl+S`).

## Memory note
MS1 spectra are kept in memory (needed for MS1-based quantification), so each imported file needs roughly 50–300 MB depending on scan count and peak density.
Remove files you no longer need (right-click → Remove from project).

## Test
```powershell
cargo test                      # unit tests + synthetic-data tests
cargo test -- --nocapture       # also shows which tests were skipped
```
Tests that use real files read them from `test-data/` (`*.raw`, optional `*.csv` target tables and `*.quant` projects). They are
skipped when the files or ThermoRawFileParser are missing. Raw files are large and are not committed (see `.gitignore`).

## Credits and licenses
- Built with **Claude Sonnet 5.5 (Medium)** (Anthropic); the **app icon was generated with Gemini** (Google).
- Thermo `.raw` reading: [ThermoRawFileParser](https://github.com/CompOmics/ThermoRawFileParser) (Apache-2.0), which uses Thermo Fisher Scientific's
  RawFileReader library under Thermo's license. Hulstaert et al., *J Proteome Res* 2020, 19, 537–542, doi:10.1021/acs.jproteome.9b00328.
- GUI: egui / eframe / egui_plot, file dialogs: rfd. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
- This project is released under the [MIT License](LICENSE).
