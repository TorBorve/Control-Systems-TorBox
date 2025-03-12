use std::env;

fn main() {
    let lib_path = format!("{}/lib", env::var("CARGO_MANIFEST_DIR").unwrap());
    println!("cargo:rustc-link-search={}", lib_path);
    println!("cargo:rustc-link-lib=gfortran");
}
