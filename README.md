# LC-MS/MS Quant

![App Screenshot](img/lcms-quant-ico.png)

Desktop app (Windows, macOS) for **targeted quantification of LC-MS/MS data** from Thermo `.raw` and mzML files.
Repository: <https://github.com/markellerhorst/lcms-quant>

You provide a target table (name, retention time or `n/a`, diagnostic MS2 product ions, optionally the MS1 precursor m/z); the program extracts the
chromatograms, integrates the peaks (boundaries adjustable by hand, values update live) and reports either **integrated areas
(arbitrary units)** or **concentrations from per-target standard curves** measured in the same project.

> Intended for method development, research and routine targeted analysis. The software provides data-integrity
> features (hash-chained audit trail, file hashes, reasons for manual changes) but is **not validated or certified** for
> GxP / 21 CFR Part 11 use. Validation of the software and of the converter version is the user's responsibility.

## Features
- Import of Thermo `.raw` (via ThermoRawFileParser) and **mzML** files (any vendor, e.g. converted with ProteoWizard msconvert); data is kept **in memory only**
- Target table: load, edit live inside the program, save under any name; optional MS1 precursor per target (`p:286.1438`)
- Quantify on the **MS2 product-ion area** or the **MS1 precursor area**; manual integration by dragging boundaries; baseline subtraction; ion-ratio checks
- **Several targets in one project**: each target has its own calibration; standards carry a concentration *per target*; replicates share a level,
  are all used by default and can be excluded one by one; levels show n, mean area and RSD. Batch assignment of selected replicate files.
- Areas-only mode, results matrix for all files, retention-time alignment (TIC cross-correlation), TIC / BPC / XIC views of one or many files
- **Quick XIC** for any m/z (MS1 or MS2) without defining a target; "Add as target" makes it permanent
- Hover read-out on chromatograms (m/z, RT, scan, precursor m/z, intensity) and a linked **MS2 spectrum panel**; the cursor snaps to the nearest peak (RT + m/z)
- Explorer-style multi-select (Ctrl / Shift), right-click to remove a file from the project (never deletes it on disk)
- Natural zoom: drag a box around a peak to zoom to it (chromatograms and MS2 spectra), Shift/middle-drag pans, right-click → *Reset zoom*;
  integration boundaries are moved by grabbing the start/end line itself
- Right-click any chromatogram or MS2 spectrum → save as **PNG**, **SVG** (spectra carry the m/z of the 10 strongest peaks) or a **two-column CSV** (RT/Intensity, m/z/Intensity); file names say what is exported (TIC, BPC, XIC, MS2 scan)
- Projects (`.quant`) that **contain the raw files** (optional), undo/redo (10 steps), zoom with Ctrl +/-/0, user name in the audit trail
- Numbers in scientific notation (`10²`), except retention times (`7.20` min) and m/z values (`286.1438`); normalised views show 0–100 %

## Build
Requirements: [Rust](https://rustup.rs) (stable); on Windows additionally the *Visual Studio Build Tools* with the
"Desktop development with C++" workload; on macOS the Xcode command-line tools.

```powershell
git clone https://github.com/markellerhorst/lcms-quant
cd lcms-quant
cargo run --release          # or: cargo build --release  -> target/release/lcms-quant(.exe)
```

### Thermo reader (only for `.raw` files)
1. Download a release of [ThermoRawFileParser](https://github.com/CompOmics/ThermoRawFileParser/releases) and unzip it to
   `tools/ThermoRawFileParser<anything>/` (searched automatically) or point to it under **Edit → Preferences**.
2. The framework-dependent builds need the **.NET 8 runtime**; the self-contained builds do not.
   Preferences → *Test* runs `ThermoRawFileParser -v`.

## Use
1. Drop a target table `.csv` and `.raw` / `.mzML` files onto the window (or use the File menu). Table format: `name,rt_min,ions` – `n/a` for an unknown RT,
   MS2 product ions separated by `;` (the first is the quantifier) and, optionally, the MS1 precursor as `p:286.1438`, e.g.
   `Morphine,n/a,286.14;201.07;165.05;p:286.1438`.
2. **Calibration setup** tab: select the replicate files of one level (Ctrl/Shift+click in the sample list), choose the target and the concentration,
   press *Assign as standards*; repeat for every level and target. Samples keep the type *Unknown*. The matrix below allows per-cell edits.
   Setting the type (Standard, Unknown, QC, Blank) of a file that is part of a multi-selection sets **all selected files** to that type.
3. Choose the signal (MS2 product ions or MS1 precursor) in the top bar. Check/adjust the integrated peaks; give a reason for manual changes.
   Replicates can be excluded from a curve with the checkbox in the calibration table.
4. Export the results (`Export CSV`), save the project (`Ctrl+S`). A project embeds the original files, so it can be opened on another PC
   (ThermoRawFileParser is needed there for `.raw`). While opening, each embedded `.raw` is copied to a temporary file for the converter and deleted again
   immediately; the data itself stays in memory. Untick *Embed raw files when saving projects* to save only references + SHA-256.

## Project file size
A `.quant` project that embeds its raw files is as large as the raw files together (often several GB). Saving runs in the background and
streams the data. Projects from older versions can still be opened.

## Memory note
MS1 spectra are kept in memory (needed for MS1-based quantification), so each imported file needs roughly 50–300 MB depending on scan count and peak density.
Remove files you no longer need (right-click → Remove from project).

## Test
```powershell
cargo test --no-fail-fast       # unit tests + synthetic-data tests (+ real files from test-data/)
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
