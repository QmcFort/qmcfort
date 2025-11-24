#!/usr/bin/env bash
set -Eeuo pipefail

basis=$1

BASH_CALL="${BASH_SOURCE[0]:-$0}"
TOOLS_PATH="$(cd "$(dirname "$BASH_CALL")" && pwd)"
BASIS_PATH="$TOOLS_PATH/../data/basis_sets"

rm -f basis_set
touch basis_set

for ((i = 2; i <= $#; i++ )); do
  cat $BASIS_PATH/$basis/${!i} >> basis_set
done
