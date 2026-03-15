use std::fs;
use std::io::Write;

use crate::simulate::SimulationResult;

// ============================================================================
// Parameter parsing
// ============================================================================

#[derive(Debug, Clone)]
pub struct SimParams {
    pub n: usize,
    pub vmode: f64,
    pub alpha_l: f64,
    pub alpha_r: f64,
    pub lambda: f64,
    pub t: f64,
    pub eta: f64,
    pub vg: f64,
    pub tau: f64,
    pub vsd_start: f64,
    pub vsd_end: f64,
    pub vsd_step: f64,
}

pub fn parse_params_json(filepath: &str) -> Option<SimParams> {
    let text = fs::read_to_string(filepath).ok()?;
    let root: serde_json::Value = serde_json::from_str(&text).ok()?;

    let params = root.get("parameters")?;
    let sweep = root.get("bias_sweep")?;

    let tau_val = match params.get("tau")? {
        serde_json::Value::String(s) if s == "Inf" => f64::INFINITY,
        serde_json::Value::Number(n) => n.as_f64().unwrap_or(f64::INFINITY),
        _ => f64::INFINITY,
    };

    Some(SimParams {
        n: params.get("N")?.as_u64()? as usize,
        vmode: params.get("vmode")?.as_f64()?,
        alpha_l: params.get("alphaL")?.as_f64()?,
        alpha_r: params.get("alphaR")?.as_f64()?,
        lambda: params.get("lambda")?.as_f64()?,
        t: params.get("T")?.as_f64()?,
        eta: params.get("eta")?.as_f64()?,
        vg: params.get("Vg")?.as_f64()?,
        tau: tau_val,
        vsd_start: sweep.get("Vsd_start")?.as_f64()?,
        vsd_end: sweep.get("Vsd_end")?.as_f64()?,
        vsd_step: sweep.get("Vsd_step")?.as_f64()?,
    })
}

// ============================================================================
// MATLAB reference loading
// ============================================================================

pub fn load_matlab_vsd(filepath: &str) -> Option<Vec<f64>> {
    let text = fs::read_to_string(filepath).ok()?;
    let root: serde_json::Value = serde_json::from_str(&text).ok()?;
    let arr = root.get("Vsd")?.as_array()?;
    Some(arr.iter().filter_map(|v| v.as_f64()).collect())
}

pub struct MatlabReference {
    pub i_tol: Vec<f64>,
    pub i_seq: Vec<f64>,
    pub i_cot: Vec<f64>,
}

pub fn load_matlab_reference(filepath: &str) -> Option<MatlabReference> {
    let text = fs::read_to_string(filepath).ok()?;
    let root: serde_json::Value = serde_json::from_str(&text).ok()?;

    let i_tol: Vec<f64> = root
        .get("I_tol")?
        .as_array()?
        .iter()
        .filter_map(|v| v.as_f64())
        .collect();
    let i_seq: Vec<f64> = root
        .get("I_seq")?
        .as_array()?
        .iter()
        .filter_map(|v| v.as_f64())
        .collect();
    let i_cot: Vec<f64> = root
        .get("I_cot")?
        .as_array()?
        .iter()
        .filter_map(|v| v.as_f64())
        .collect();

    Some(MatlabReference { i_tol, i_seq, i_cot })
}

// ============================================================================
// Output writing
// ============================================================================

fn write_json_array(f: &mut fs::File, name: &str, arr: &[f64], comma: bool) -> std::io::Result<()> {
    write!(f, "  \"{}\": [\n", name)?;
    for (i, &v) in arr.iter().enumerate() {
        if i + 1 < arr.len() {
            writeln!(f, "    {},", v)?;
        } else {
            writeln!(f, "    {}", v)?;
        }
    }
    if comma {
        writeln!(f, "  ],")?;
    } else {
        writeln!(f, "  ]")?;
    }
    Ok(())
}

pub fn write_results_json(
    filepath: &str,
    spec: &str,
    params: &SimParams,
    result: &SimulationResult,
    wall_time: f64,
) -> std::io::Result<()> {
    let mut f = fs::File::create(filepath)?;

    // Get timestamp
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    // Simple timestamp formatting (UTC)
    let secs_per_day = 86400u64;
    let secs_per_hour = 3600u64;
    let secs_per_min = 60u64;

    // Days since epoch to Y-M-D
    let mut days = now / secs_per_day;
    let time_of_day = now % secs_per_day;
    let hour = time_of_day / secs_per_hour;
    let minute = (time_of_day % secs_per_hour) / secs_per_min;
    let second = time_of_day % secs_per_min;

    // Simple year/month/day calculation
    let mut year = 1970i32;
    loop {
        let days_in_year = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
            366
        } else {
            365
        };
        if days < days_in_year {
            break;
        }
        days -= days_in_year;
        year += 1;
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days_in_months: [u64; 12] = [
        31,
        if leap { 29 } else { 28 },
        31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
    ];
    let mut month = 1u32;
    for &dim in &days_in_months {
        if days < dim {
            break;
        }
        days -= dim;
        month += 1;
    }
    let day = days + 1;

    let timestamp = format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
        year, month, day, hour, minute, second
    );

    let rust_version = env!("CARGO_PKG_VERSION");

    writeln!(f, "{{")?;
    writeln!(f, "  \"spec\": \"{}\",", spec)?;
    writeln!(f, "  \"language\": \"Rust\",")?;
    writeln!(f, "  \"version\": \"rustc (franck-condon {})\",", rust_version)?;
    writeln!(f, "  \"wall_time_seconds\": {},", wall_time)?;
    writeln!(f, "  \"parameters\": {{")?;
    writeln!(f, "    \"N\": {},", params.n)?;
    writeln!(f, "    \"vmode\": {},", params.vmode)?;
    writeln!(f, "    \"alphaL\": {},", params.alpha_l)?;
    writeln!(f, "    \"alphaR\": {},", params.alpha_r)?;
    writeln!(f, "    \"lambda\": {},", params.lambda)?;
    writeln!(f, "    \"T\": {},", params.t)?;
    writeln!(f, "    \"eta\": {},", params.eta)?;
    writeln!(f, "    \"Vg\": {},", params.vg)?;
    if params.tau.is_infinite() {
        writeln!(f, "    \"tau\": \"Inf\"")?;
    } else {
        writeln!(f, "    \"tau\": {}", params.tau)?;
    }
    writeln!(f, "  }},")?;
    writeln!(f, "  \"bias_sweep\": {{")?;
    writeln!(f, "    \"Vsd_start\": {},", params.vsd_start)?;
    writeln!(f, "    \"Vsd_end\": {},", params.vsd_end)?;
    writeln!(f, "    \"Vsd_step\": {}", params.vsd_step)?;
    writeln!(f, "  }},")?;

    write_json_array(&mut f, "Vsd", &result.vsd, true)?;
    write_json_array(&mut f, "I_tol", &result.i_tol, true)?;
    write_json_array(&mut f, "I_seq", &result.i_seq, true)?;
    write_json_array(&mut f, "I_cot", &result.i_cot, true)?;
    writeln!(f, "  \"timestamp\": \"{}\"", timestamp)?;
    writeln!(f, "}}")?;

    Ok(())
}

pub fn write_results_csv(filepath: &str, result: &SimulationResult) -> std::io::Result<()> {
    let mut f = fs::File::create(filepath)?;
    writeln!(f, "Vsd_V,I_tol_A,I_seq_A,I_cot_A")?;
    for i in 0..result.vsd.len() {
        writeln!(
            f,
            "{:.6e},{:.6e},{:.6e},{:.6e}",
            result.vsd[i], result.i_tol[i], result.i_seq[i], result.i_cot[i]
        )?;
    }
    Ok(())
}
