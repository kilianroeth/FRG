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

printFRG() {
    local C1="\033[38;2;50;255;120m"  # Emerald Green
    local C2="\033[38;2;20;240;160m"  # Mint Green
    local C3="\033[38;2;0;225;200m"   # Bright Cyan
    local C4="\033[38;2;0;200;230m"   # Ocean Cyan
    local C5="\033[38;2;0;160;255m"   # Sky Blue
    local C6="\033[38;2;30;110;255m"  # Royal Blue
    local C7="\033[38;2;60;60;255m"   # Deep Electric Blue
    local RESET="\033[0m"

    printf "=================================\n"
    printf "\n"
    printf "${C1}   ███████╗██████╗  ██████╗  ${RESET}\n"
    printf "${C2}   ██╔════╝██╔══██╗██╔════╝  ${RESET}\n"
    printf "${C3}   █████╗  ██████╔╝██║  ███╗ ${RESET}\n"
    printf "${C4}   ██╔══╝  ██╔══██╗██║   ██║ ${RESET}\n"
    printf "${C5}   ██║     ██║  ██║╚██████╔╝ ${RESET}\n"
    printf "${C6}   ╚═╝     ╚═╝  ╚═╝ ╚═════╝  ${RESET}\n"
    printf "\n"
    printf "==================================\n"
}

for arg in "$@"; do
    case $arg in
        --run|-r)      RUN=true ;;
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

printFRG

if $RUN; then
    echo ""

    # Determine 2/3 of available logical cores (rounded down, minimum 1)
    TOTAL_CORES=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)
    OMP_THREADS=$(( TOTAL_CORES * 5 / 6 ))
    [[ "$OMP_THREADS" -lt 1 ]] && OMP_THREADS=1

    info "Running $EXECUTABLE with OMP_NUM_THREADS=$OMP_THREADS (of $TOTAL_CORES cores)..."
    echo "----------------------------------"
    export OMP_NUM_THREADS=$OMP_THREADS
    "$EXECUTABLE"
fi