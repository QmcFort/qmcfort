#!/usr/bin/env bash
set -Eeuo pipefail

# ==============================================================================
# QmcFort regression tests driver:
#   - Test selection with globs (e.g., "*afqm*c" "*hf*" "benzene_afqmc_vdz")
#   - Run tests in fresh working copies under results/<case>
#   - Updates tests by calling per-test update_test()
# ==============================================================================


# Locate relevant paths relative to this file ==================================
REG_SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REG_ROOT="$(cd "${REG_SRC_DIR}/.." && pwd)"         # qmcfort/tests/regression/
export REG_ROOT 


# Defaults (can be set in config.sh) ===========================================
RESULTS_DIR_DEFAULT="${RESULTS_DIR_DEFAULT:-${REG_ROOT}/results}"
CASES_DEFAULT="${CASES_DEFAULT:-}"


# Sync two directories =========================================================
function sync_regression_test_dir() {
  local src="$1" 
  local dest="$2"

  rm -rf "$dest"
  mkdir -p "$dest"

  rsync -a --delete "$src/" "$dest/"
}


# Test selection logic =========================================================
# Usage:
#   select_tests "$CASES"              # patterns like "*afqmc*", "benzene_vdz"
#   select_tests ""                    # select all
# Populates global array: SELECTED
# ==============================================================================
function select_tests() {
  local cases_arg="${1:-}"
  local tests_dir="${REG_ROOT}/cases"

  # 1) discover all test folders (names only)
  mapfile -t ALL_CASES < <(find "${tests_dir}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

  # 2) if no filters given -> select all
  if [[ -z "$cases_arg" ]]; then
    SELECTED=("${ALL_CASES[@]}")
    return 0
  fi

  # 3) split patterns (space-separated)
  local -a PATS
  read -r -a PATS <<<"$cases_arg"

  # 4) match with shell globs; deduplicate via assoc array
  declare -A seen=()
  local c p
  SELECTED=()
  for c in "${ALL_CASES[@]}"; do
    for p in "${PATS[@]}"; do
      # glob match (supports * ? […])
      if [[ "$c" == $p ]]; then
        if [[ -z "${seen[$c]:-}" ]]; then
          SELECTED+=("$c")
          seen[$c]=1
        fi
        break
      fi
    done
  done

  if [[ ${#SELECTED[@]} -eq 0 ]]; then
    echo "[WARN] No tests matched patterns: ${PATS[*]}"
    return 1
  fi
  return 0
}


# Wrapper around cases/<case>/test.sh:run_test() ===============================
#   - create fresh working copy in out_dir/<case>
#   - execute cases/<case>/test.sh:run_test() inside the working dir.
# ==============================================================================
run_regression_test() {
  local out_dir=$1
  local case=$2
  local test_counter=$3

  local src="${REG_ROOT}/cases/${case}"
  local dest="${out_dir}/${case}"

  sync_regression_test_dir "$src" "$dest"

  (
    cd "$dest"

    test_intro $case $test_counter

    # Execute the test rules
    source "./test.sh"
    run_test
  )
}


# Run selected regression tests ================================================
# Usage: 
#   run_regression_tests build/bin/qmcfort tests/regression/results "*afqmc*"
# ==============================================================================
function run_regression_tests() {
  local bin=$1
  local out_dir="${2:-$RESULTS_DIR_DEFAULT}"
  local cases="${3:-$CASES_DEFAULT}"

  export BIN=$bin

  [[ -n "$BIN" && -x "$BIN" ]] || { echo "[ERR] --bin must be an executable"; return 3; }

  mkdir -p "$out_dir"

  if ! select_tests "$cases"; then
    echo "[INFO] Nothing to run."
    return 0
  fi

  # Loop over cases and their indices
  echo "[INFO] Running ${#SELECTED[@]} test(s)"
  for i in "${!SELECTED[@]}"
  do
    echo "[INFO]    $((i + 1)).    ${SELECTED[$i]}"
  done

  local test_counter=0

  python3 $REG_ROOT/src/timestamp.py

  for case in "${SELECTED[@]}"
  do
    test_counter=$(($test_counter + 1))
    run_regression_test "$out_dir" "$case" "$test_counter"
  done

  python3 $REG_ROOT/src/check_report.py $REG_ROOT/report.log
}


# Update ref files of selected regression tests ================================
#   - Calls each test's update_test() inside WORK, 
#   - sync out_dir/case/ref_data -> tests/case/ref_data
# Usage: 
#   update_regression_tests "reg_tests/results" "benzene_afqmc_vdz h2o_hf_vqz"
# ==============================================================================
function update_regression_tests() {
  local out_dir="${1:-$RESULTS_DIR_DEFAULT}"
  local cases="${2:-$CASES_DEFAULT}"

  if ! select_tests "$cases"; then
    echo "[INFO] Nothing to update."
    return 0
  fi

  for case in "${SELECTED[@]}"
  do
    local dest="$out_dir/$case"
    local src="${REG_ROOT}/cases/$case"

    (
      cd "$dest"

      # Load original test rules from source tree
      source "$src/test.sh"

      if ! declare -F update_test >/dev/null; then
        echo "[WARN] $case: update_test() not defined in test.sh"
        exit 11
      fi

      # Let the test implement its own sync policy
      update_test "$dest" "$src"
    )
  done
}


function test_intro() {
    local test_name=$1
    local test_counter=$2

    echo "################################################################################"
    echo ""
    echo "Regression test number $test_counter:  $test_name"
    echo ""
    echo "################################################################################"
}