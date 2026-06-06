#!/bin/bash

set -euo pipefail

CYAN="\e[36m"
GREEN="\e[32m"
RED="\e[31m"
RESET="\e[0m"

info() { echo -e "$CYAN[INFO]$RESET $*"; }
error() { echo -e "$RED[ERROR]$RESET $*" >&2; exit 1; }
success() { echo -e "$GREEN[OK]$RESET $*"; }

BUILD_TYPE="${BUILD_TYPE:-Debug}"
BUILD_DIR="build"
RUN=false
CLEAN=false

usage() {
    echo "Usage: ./build.sh [options]"
    echo ""
    echo "Options:"
    echo "--run        Build then run the exectuable"
    echo "--release    Build in Release mode (default: Debug)"
    echo "--clean      Remove build directory before building"
    echo "-h, --help   Show this help message"
    exit 0
}

for arg in "$@"; do
    case $arg in
        --run)      RUN=true ;;
        --release)  BUILD_TYPE="Release" ;;
        --clean)    CLEAN=true ;;
        -h|--help)  usage;;
        *) error "Unknown argument: $arg" ;;
    esac
done

command -v cmake &>/dev/null || error "cmake not found. Please install CMake."

if $CLEAN; then
    info "Cleaing build directory..."
    rm -rf "$BUILD_DIR"
    success "Cleaned."
fi

mkdir -p "$BUILD_DIR"
info "Configuring ($BUILD_TYPE)..."
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -G Ninja

info "Building..."
cmake --build "$BUILD_DIR" --parallel "$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"

# Find the build executable
EXECUTABLE=$(find "$BUILD_DIR" -maxdepth 2 -type f -perm /111 ! -name "*.cmake" ! -name "Makefile" | head -n 1)
[[ -n "$EXECUTABLE" ]] || error "Could not find build exectuable."

success "Build $EXECUTABLE"

if $RUN; then
    echo ""
    info "Running $EXECUTABLE..."
    echo "----------------------------------"
    "$EXECUTABLE"
fi
