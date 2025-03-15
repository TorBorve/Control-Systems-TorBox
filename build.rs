use make_cmd::make;
use std::{env, path::PathBuf};

fn main() {
    let slicot_dir = PathBuf::from("SLICOT-Reference/");

    let make_file = if cfg!(target_os = "windows") {
        "makefile"
    } else {
        "makefile_Unix"
    };

    let mut make_args: Vec<String> = vec![
        "lib".to_string(),
        "-f".to_string(),
        make_file.to_string(),
        "-j8".to_string(),
    ];
    if cfg!(target_os = "windows") {
        make_args.push("HOME=.".to_string());
    } else {
        make_args.push("SLICOTLIB=../libslicot.a".to_string());
        make_args.push("LPKAUXLIB=../liblpkaux.a".to_string());
        if let Ok(fc) = env::var("FC") {
            make_args.push(format!("FORTRAN={}", fc).to_string());
        }
    }

    let status = make()
        .current_dir(&slicot_dir)
        .args(make_args)
        .status()
        .expect("Failed to build make");

    if !status.success() {
        panic!("Failed to build SLICOT")
    }

    println!("cargo:rustc-link-search={}", slicot_dir.to_str().unwrap());
    println!("cargo:rustc-link-lib=static=slicot");
    println!("cargo:rustc-link-lib=static=lpkaux");
    println!("cargo:rustc-link-lib=gfortran");
}
