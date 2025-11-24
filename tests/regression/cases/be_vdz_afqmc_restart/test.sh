#!/usr/bin/env bash
set -Eeuo pipefail

# Resolve directory names (TEST_NAME is just a name of the test directory)
REG_TEST_PATH="${BASH_SOURCE[0]:-$0}"
REG_TEST_DIR="$(cd "$(dirname "$REG_TEST_PATH")" && pwd)"
REG_TEST_ROOT="$(cd "$REG_TEST_DIR/../.." && pwd)"
TEST_NAME="$(basename "$REG_TEST_DIR")"

source $REG_TEST_ROOT/src/regression_test.sh
source $REG_TEST_ROOT/src/extractors.sh

function run_test() {
  cp qmcfort_in_hf qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_hf
  cp qmcfort_log qmcfort_log_hf
  cp timing timing_hf

  cp qmcfort_in_afqmc1 qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_afqmc1
  cp qmcfort_log qmcfort_log_afqmc1
  cp timing timing_afqmc1

  cp qmcfort_in_afqmc2 qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_afqmc2
  cp qmcfort_log qmcfort_log_afqmc2
  cp timing timing_afqmc2

  #extract values for comparison
  rm -f values.txt
  touch values.txt
  afqmc_search "qmcfort_out_afqmc2" "values.txt"

  #compare values.txt and ref_values.txt
  python3 $REG_TEST_ROOT/src/test_report.py $TEST_NAME
}

function clean_test() {
  qmcfort_clean
}

function update_test() {
  local result="$1"       # reg_tests/results/<test>
  local source="$2"       # reg_tests/tests/<test>

  rsync -a $result/values.txt $source/ref_values.txt
  rsync -a $result/qmcfort_out_hf $source/ref_qmcfort_out_hf.txt
  rsync -a $result/qmcfort_out_afqmc1 $source/ref_qmcfort_out_afqmc1.txt
  rsync -a $result/qmcfort_out_afqmc2 $source/ref_qmcfort_out_afqmc2.txt
  rsync -a $result/timing_hf $source/ref_timing_hf.txt
  rsync -a $result/timing_afqmc1 $source/ref_timing_afqmc1.txt
  rsync -a $result/timing_afqmc2 $source/ref_timing_afqmc2.txt
}
