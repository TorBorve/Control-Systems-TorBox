#!/usr/bin/env bash

# This script is used to test the CI/CD pipeline
set -e  # Exit immediately if a command fails

run_command() {
    echo "Running: $1"
    if ! eval "$1"; then
        echo "❌ Error: Command '$1' failed."
        exit 1
    fi
}

run_command "cargo fmt --all -- --check"
run_command "cargo clippy --all -- -D warnings"
run_command "cargo test --all"
run_command "cargo check"
run_command "cargo llvm-cov --lcov --output-path target/llvm-cov-target/lcov.info"

echo "✅ All checks passed successfully!"