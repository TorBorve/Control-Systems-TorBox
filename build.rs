use make_cmd::make;
use std::{
    env,
    path::{Path, PathBuf},
};

fn main() {
    let slicot_dir = PathBuf::from("SLICOT-Reference/");
    let slicot_build_dir = Path::new(&env::var("OUT_DIR").unwrap())
        .join("slicot_build")
        .to_str()
        .unwrap()
        .to_string();
    let make_args: Vec<String> = vec![
        "-j8".to_string(),
        format!("BUILD_DIR={slicot_build_dir}").to_string(),
    ];

    let status = make()
        .current_dir(&slicot_dir)
        .args(make_args)
        .status()
        .expect("Failed to build make");

    if !status.success() {
        if cfg!(target_os = "windows") {
            println!(
                "cargo:warning=Control Systems Torbox has not been tested on windows. You are likely to have issues with compiling SLICOT."
            )
        }
        println!("cargo:warning=Failed to compile SLICOT library");
        panic!("Failed to build SLICOT")
    }

    println!("cargo:rustc-link-search={slicot_build_dir}");
    println!("cargo:rustc-link-lib=static=slicot");
    println!("cargo:rustc-link-lib=static=lpkaux");
    println!("cargo:rustc-link-lib=gfortran");
}
