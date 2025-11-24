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
  cp qmcfort_in_ints qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_ints
  cp timing timing_ints

  cp qmcfort_in_chol qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_chol
  cp timing timing_chol

  cp qmcfort_in_svd qmcfort_in
  qmcfort_run 1 4
  cp qmcfort_out qmcfort_out_svd
  cp timing timing_svd

  #extract values for comparison
  rm -f values.txt
  touch values.txt
  hf_search "qmcfort_out_svd" "values.txt"

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
  rsync -a $result/qmcfort_out_ints $source/ref_qmcfort_out_ints.txt
  rsync -a $result/qmcfort_out_chol $source/ref_qmcfort_out_chol.txt
  rsync -a $result/qmcfort_out_svd $source/ref_qmcfort_out_svd.txt
  rsync -a $result/timing_ints $source/ref_timing_ints.txt
  rsync -a $result/timing_chol $source/ref_timing_chol.txt
  rsync -a $result/timing_svd $source/ref_timing_svd.txt
}