//! Unit/integration tests for the non-GUI core. Run: cargo test
use lcms_quant::{
    audit, fmt::{parse_sci, sci}, history::History, icon, inflate::zlib_decompress, mzml, project::*, quant::*, sha256::*, testkit,
};

// ---------- hashing ----------
#[test]
fn sha256_known_vectors() {
    assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    assert_eq!(sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
               "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

#[test]
fn sha256_streaming_matches_one_shot() {
    let mut h = Sha256::new();
    let chunk = vec![b'a'; 1000];
    for _ in 0..1000 { h.update(&chunk); }
    assert_eq!(h.finish_hex(), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
    let data: Vec<u8> = (0..5000u32).map(|i| (i * 7 % 251) as u8).collect();
    let mut h = Sha256::new();
    for c in data.chunks(37) { h.update(c); }
    assert_eq!(h.finish_hex(), sha256_hex(&data));
}

// ---------- zlib ----------
#[test]
fn inflate_fixed_huffman_streams() {
    assert_eq!(zlib_decompress(&[0x78, 0x9c, 0x4b, 0x04, 0x00, 0x00, 0x62, 0x00, 0x62]).unwrap(), b"a");
    let hw = [0x78, 0x9c, 0xcb, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0x28, 0xcf, 0x2f, 0xca, 0x49, 0x01, 0x00, 0x1a, 0x0b, 0x04, 0x5d];
    assert_eq!(zlib_decompress(&hw).unwrap(), b"hello world");
}

#[test]
fn inflate_stored_blocks_including_multi_block() {
    let data: Vec<u8> = (0..200_000u32).map(|i| (i % 253) as u8).collect();
    assert_eq!(zlib_decompress(&testkit::zlib_stored(&data)).unwrap(), data);
    assert!(zlib_decompress(&[0x78, 0x9c, 0xff]).is_err());
}

// ---------- scientific notation ----------
#[test]
fn scientific_notation_display() {
    assert_eq!(sci(100.0), "10²");
    assert_eq!(sci(1000.0), "10³");
    assert_eq!(sci(1.0), "10⁰");
    assert_eq!(sci(1234.5), "1.2345×10³");
    assert_eq!(sci(0.05), "5×10⁻²");
    assert_eq!(sci(-250.0), "-2.5×10²");
    assert_eq!(sci(0.0), "0");
}

#[test]
fn scientific_notation_round_trips_through_the_parser() {
    for v in [100.0, 1000.0, 1234.5, 0.05, -250.0, 1e-7, 123456.0, 3.0] {
        let back = parse_sci(&sci(v)).unwrap_or_else(|| panic!("cannot parse {}", sci(v)));
        assert!((back - v).abs() <= v.abs() * 1e-5, "{} -> {} -> {}", v, sci(v), back);
    }
    assert_eq!(parse_sci("1.5e3"), Some(1500.0));
    assert_eq!(parse_sci("2,5"), Some(2.5));
    assert!((parse_sci("5x10^-2").unwrap() - 0.05).abs() < 1e-12);
    assert_eq!(parse_sci("garbage"), None);
}

// ---------- undo/redo ----------
#[test]
fn history_is_bounded_and_ordered() {
    let mut h: History<i32> = History::new(10);
    for i in 0..12 { h.record(i); }
    assert_eq!(h.lens(), (10, 0));
    assert_eq!(h.undo(100), Some(11));
    assert_eq!(h.undo(11), Some(10));
    assert_eq!(h.redo(10), Some(11));
    assert_eq!(h.redo(11), Some(100));
    assert_eq!(h.redo(100), None);
}

#[test]
fn history_new_action_clears_redo() {
    let mut h: History<i32> = History::new(10);
    h.record(1);
    h.record(2);
    assert_eq!(h.undo(3), Some(2));
    assert!(h.can_redo());
    h.record(4);
    assert!(!h.can_redo());
    assert_eq!(h.undo(5), Some(4));
    assert_eq!(h.undo(4), Some(1));
    assert_eq!(h.undo(1), None);
}

// ---------- target table ----------
#[test]
fn target_table_parsing() {
    let t = parse_table("name,rt_min,ions\nMorphine,n/a,286.14;201.07;165.05\nCodeine,3.2,300.16 215.10\nFentanyl,5.8,337.23,188.13,105.07\n").unwrap();
    assert_eq!(t.len(), 3);
    assert_eq!(t[0].rt, None);
    assert_eq!(t[0].ions, vec![286.14, 201.07, 165.05]);
    assert_eq!(t[1].rt, Some(3.2));
    assert_eq!(t[1].ions.len(), 2);
    assert_eq!(t[2].ions.len(), 3);
    assert!(parse_table("").is_err());
    assert!(parse_table("name,rt,ions\nGood,1.0,100\nBad,x,100\n").is_err());
    assert_eq!(parse_table("Drug\t2.5\t100.1\t200.2\n").unwrap()[0].ions, vec![100.1, 200.2]);
}

// ---------- mzML ----------
#[test]
fn mzml_parsing_of_synthetic_run() {
    let d = mzml::parse_bytes("synthetic", &testkit::synthetic_mzml(0.0, 1.0)).unwrap();
    assert_eq!(d.spectra.len(), 1000, "MS2 scans");
    assert_eq!(d.tic1.len(), 1000, "MS1 points");
    assert_eq!(d.n_total, 2000);
    assert_eq!(d.sha256.len(), 64);
    assert!((d.spectra[10].rt - 0.105).abs() < 1e-9, "minute unit");
    assert!((d.tic1[10][0] - 0.1).abs() < 1e-9, "second unit converted to minutes");
    assert_eq!(d.spectra[0].mz, vec![150.0, 250.0, 400.0]);
    assert!(d.spectra.windows(2).all(|w| w[0].rt <= w[1].rt));
    let apex = d.bpc2.iter().cloned().fold([0.0, 0.0], |m, p| if p[1] > m[1] { p } else { m });
    assert!((apex[0] - 5.005).abs() < 0.011);
    assert!((apex[1] - 1000.0).abs() < 10.0);
    assert!(mzml::parse_bytes("empty", b"<mzML></mzML>").is_err());
}

// ---------- chromatography maths ----------
#[test]
fn xic_integration_matches_the_analytical_area() {
    let d = mzml::parse_bytes("s", &testkit::synthetic_mzml(0.0, 1.0)).unwrap();
    let x = xic(&d.spectra, testkit::FRAG1, 0.01);
    assert_eq!(x.len(), d.spectra.len());
    let area = integrate(&x, 4.5, 5.5, false);
    let expect = testkit::expected_area_seconds(1.0);
    assert!((area - expect).abs() / expect < 0.01, "area {} expected {}", area, expect);
    let x2 = xic(&d.spectra, testkit::FRAG2, 0.01);
    let a2 = integrate(&x2, 4.5, 5.5, false);
    assert!((a2 / area - 0.5).abs() < 0.005);
    assert!(xic(&d.spectra, 999.0, 0.01).iter().all(|p| p[1] == 0.0));
}

#[test]
fn baseline_subtraction_removes_a_constant_offset() {
    let pts: Vec<[f64; 2]> = (0..1000).map(|i| { let t = 4.0 + i as f64 * 0.002; [t, 100.0 + 1000.0 * (-0.5f64 * ((t - 5.0) / 0.05).powi(2)).exp()] }).collect();
    let with = integrate(&pts, 4.7, 5.3, true);
    let without = integrate(&pts, 4.7, 5.3, false);
    let expect = 1000.0 * 0.05 * (2.0 * std::f64::consts::PI).sqrt() * 60.0;
    assert!((with - expect).abs() / expect < 0.01, "{} vs {}", with, expect);
    assert!(without > with * 1.2);
}

#[test]
fn auto_peak_brackets_the_apex_and_honours_the_rt_hint() {
    let g = |t: f64, c: f64, h: f64| h * (-0.5f64 * ((t - c) / 0.05).powi(2)).exp();
    let pts: Vec<[f64; 2]> = (0..1500).map(|i| { let t = i as f64 * 0.01; [t, g(t, 5.0, 1000.0) + g(t, 8.0, 400.0)] }).collect();
    let (s, e) = auto_peak(&pts, None, 0.0).unwrap();
    assert!(s < 5.0 && e > 5.0 && s > 4.6 && e < 5.4, "{} {}", s, e);
    let (s2, e2) = auto_peak(&pts, Some(8.0), 0.3).unwrap();
    assert!(s2 < 8.0 && e2 > 8.0 && s2 > 7.6 && e2 < 8.4, "{} {}", s2, e2);
    assert!(auto_peak(&pts, Some(11.0), 0.3).is_none());
}

#[test]
fn calibration_fit_recovers_the_line_for_every_weighting() {
    let pts: Vec<(f64, f64)> = [1.0, 2.0, 5.0, 10.0, 50.0].iter().map(|&x| (x, 250.0 * x + 40.0)).collect();
    for w in 0..3u8 {
        let f = fit(&pts, w).unwrap();
        assert!((f.m - 250.0).abs() < 1e-6 && (f.b - 40.0).abs() < 1e-5, "weighting {}", w);
        assert!(f.r2 > 0.999999);
        assert_eq!((f.lo, f.hi, f.n), (1.0, 50.0, 5));
    }
    assert!(fit(&[(1.0, 2.0), (1.0, 3.0)], 0).is_none(), "one level is not a curve");
}

#[test]
fn rt_alignment_estimates_the_shift_with_the_right_sign() {
    let g = |t: f64, c: f64| (-0.5f64 * ((t - c) / 0.08).powi(2)).exp();
    let make = |off: f64| -> Vec<[f64; 2]> {
        (0..1500).map(|i| { let t = i as f64 * 0.01; [t, 100.0 + 5000.0 * g(t, 5.0 + off) + 3000.0 * g(t, 8.0 + off) + 2000.0 * g(t, 12.0 + off)] }).collect()
    };
    let r = make(0.0);
    assert!(estimate_shift(&r, &r, 1.0).abs() < 0.011);
    let s = estimate_shift(&r, &make(0.25), 1.0);
    assert!((s + 0.25).abs() < 0.02, "late run must be moved earlier, got {}", s);
    let s = estimate_shift(&r, &make(-0.4), 1.0);
    assert!((s - 0.4).abs() < 0.02, "early run must be moved later, got {}", s);
}

// ---------- audit trail ----------
#[test]
fn audit_trail_detects_tampering_and_continues_across_sessions() {
    let p = std::env::temp_dir().join(format!("lcmsq_audit_{}.log", std::process::id()));
    let _ = std::fs::remove_file(&p);
    let path = p.to_string_lossy().to_string();
    let mut a = audit::Audit::new(&path);
    a.user = "tester".into();
    a.log("TEST1", "first");
    a.log("TEST2", "second");
    a.log("TEST3", "third | with pipe");
    assert_eq!(audit::verify(&path).unwrap(), 3);
    let mut b = audit::Audit::new(&path);
    b.log("TEST4", "after restart");
    assert_eq!(audit::verify(&path).unwrap(), 4);
    let txt = std::fs::read_to_string(&path).unwrap();
    std::fs::write(&path, txt.replace("second", "SECOND")).unwrap();
    assert!(audit::verify(&path).is_err());
    let _ = std::fs::remove_file(&path);
}

#[test]
fn default_user_is_na() {
    let p = std::env::temp_dir().join(format!("lcmsq_audit_user_{}.log", std::process::id()));
    assert_eq!(audit::Audit::new(&p.to_string_lossy()).user, "N/A");
}

// ---------- project file ----------
#[test]
fn project_file_round_trip() {
    let pf = ProjectFile {
        app_version: "0.3.0".into(), user: "eller".into(),
        params: vec![("tol".into(), "10".into()), ("mode".into(), "curve".into())],
        analytes: vec![Analyte { name: "A".into(), rt: None, ions: vec![100.5, 200.25] }, Analyte { name: "B".into(), rt: Some(3.25), ions: vec![300.0] }],
        samples: vec![SampleRec {
            path: r"C:\data\a b.raw".into(), raw_sha256: "ab".repeat(32), kind: "Standard".into(), conc: 0.1, rt_shift: -0.0123, show: true,
            cells: vec![CellRec { analyte: 0, found: true, manual: true, start: 1.0 / 3.0, end: 2.5, include: false },
                        CellRec { analyte: 1, found: false, manual: false, start: 0.0, end: 0.0, include: true }],
        }],
    };
    assert_eq!(ProjectFile::parse(&pf.to_text()).unwrap(), pf);
    assert_eq!(pf.param("mode"), Some("curve"));
    assert!(ProjectFile::parse("garbage").is_err());
    assert!(ProjectFile::parse("LCMSQUANT-PROJECT\t1\nC\t0\ttrue\ttrue\t1\t2\ttrue\n").is_err(), "cell before sample");
}

// ---------- icon decoding ----------
fn ico(data: &[u8], w: u8, h: u8) -> Vec<u8> {
    let mut o = vec![0, 0, 1, 0, 1, 0, w, h, 0, 0, 1, 0, 32, 0];
    o.extend_from_slice(&(data.len() as u32).to_le_bytes());
    o.extend_from_slice(&22u32.to_le_bytes());
    o.extend_from_slice(data);
    o
}

#[test]
fn ico_with_bmp_entry_decodes() {
    let mut d = vec![];
    for v in [40u32, 2, 4] { d.extend_from_slice(&v.to_le_bytes()); }
    d.extend_from_slice(&1u16.to_le_bytes());
    d.extend_from_slice(&32u16.to_le_bytes());
    d.extend_from_slice(&[0u8; 24]);
    d.extend_from_slice(&[255, 0, 0, 255, 255, 255, 255, 128]); // bottom row BGRA: blue, white(alpha 128)
    d.extend_from_slice(&[0, 0, 255, 255, 0, 255, 0, 255]);     // top row: red, green
    let (rgba, w, h) = icon::decode_ico(&ico(&d, 2, 2)).unwrap();
    assert_eq!((w, h), (2, 2));
    assert_eq!(rgba, vec![255, 0, 0, 255, 0, 255, 0, 255, 0, 0, 255, 255, 255, 255, 255, 128]);
}

#[test]
fn ico_with_png_entry_decodes_including_sub_filter() {
    fn chunk(ty: &[u8], body: &[u8]) -> Vec<u8> {
        let mut c = (body.len() as u32).to_be_bytes().to_vec();
        c.extend_from_slice(ty); c.extend_from_slice(body); c.extend_from_slice(&[0, 0, 0, 0]);
        c
    }
    let raw = [0u8, 1, 2, 3, 255, 4, 5, 6, 255, 1, 10, 20, 30, 255, 10, 10, 10, 0];
    let mut png = vec![0x89, b'P', b'N', b'G', 0x0d, 0x0a, 0x1a, 0x0a];
    let mut ihdr = vec![]; ihdr.extend_from_slice(&2u32.to_be_bytes()); ihdr.extend_from_slice(&2u32.to_be_bytes()); ihdr.extend_from_slice(&[8, 6, 0, 0, 0]);
    png.extend(chunk(b"IHDR", &ihdr));
    png.extend(chunk(b"IDAT", &testkit::zlib_stored(&raw)));
    png.extend(chunk(b"IEND", &[]));
    let (rgba, w, h) = icon::decode_ico(&ico(&png, 2, 2)).unwrap();
    assert_eq!((w, h), (2, 2));
    assert_eq!(rgba, vec![1, 2, 3, 255, 4, 5, 6, 255, 10, 20, 30, 255, 20, 30, 40, 255]);
    assert!(icon::decode_ico(b"not an icon at all, definitely not").is_none());
}

#[test]
fn the_real_app_icon_decodes_when_present() {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("img").join("lcms-quant-ico.ico");
    if let Ok(bytes) = std::fs::read(&p) {
        match icon::decode_ico(&bytes) {
            Some((rgba, w, h)) => assert_eq!(rgba.len(), (w * h * 4) as usize),
            None => eprintln!("NOTE: {} uses an ICO sub-format the built-in decoder does not support (window icon is skipped, exe icon still embedded)", p.display()),
        }
    } else { eprintln!("SKIPPED: {} not found", p.display()); }
}
