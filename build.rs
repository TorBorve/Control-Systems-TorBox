use make_cmd::make;
use std::path::PathBuf;

fn main() {
    let slicot_dir = PathBuf::from("SLICOT-Reference/");

    let make_file = if cfg!(target_os = "windows") {
        "makefile"
    } else {
        "makefile_Unix"
    };

    let mut make_args = vec!["lib", "-f", make_file, "-j8"];
    if !cfg!(target_os = "windows") {
        make_args.push("LPKAUXLIB=../liblpkaux.a");
        make_args.push("SLICOTLIB=../libslicot.a");
    } else {
        make_args.push("HOME=.");
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
