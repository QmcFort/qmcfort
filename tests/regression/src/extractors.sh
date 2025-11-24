#!/usr/bin/env bash
set -Eeuo pipefail

#
# Extract important values from output files for each of the QmcFort methods
#

function hf_search() {
    local OUT_FILE=$1
    local VAL_FILE=$2

    awk '/hf_enuc /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/hf_e1   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/hf_eh   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/hf_ex   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/hf_e    /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/hf_etot /{print $3}' $OUT_FILE  >> $VAL_FILE
}

function mp2_search() {
    local OUT_FILE=$1
    local VAL_FILE=$2

    awk '/mp2_enuc /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/mp2_e1   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/mp2_eh   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/mp2_ex   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/mp2_e    /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/mp2_etot /{print $3}' $OUT_FILE  >> $VAL_FILE
}

function afqmc_search() {
    local OUT_FILE=$1
    local VAL_FILE=$2

    awk '/trial_enuc /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/trial_e1   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/trial_eh   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/trial_ex   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/trial_e    /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/trial_etot /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_enuc /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_e1   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_e1   /{print $7}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_eh   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_eh   /{print $7}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_ex   /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_ex   /{print $7}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_e    /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_e    /{print $7}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_etot /{print $3}' $OUT_FILE  >> $VAL_FILE
    awk '/afqmc_etot /{print $7}' $OUT_FILE  >> $VAL_FILE
}

function noci_afqmc_search() {
    local OUT_FILE=$1
    local VAL_FILE=$2

    awk '/noci:f:/{print $4}' $OUT_FILE >> $VAL_FILE
    awk '/noci:f:/{print $5}' $OUT_FILE >> $VAL_FILE
    awk '/noci:f:/{print $6}' $OUT_FILE >> $VAL_FILE
}

function safqmc_search() {
    local OUT_FILE=$1
    local VAL_FILE=$2

    awk '/safqmc:/{print $3}' $OUT_FILE >> $VAL_FILE
    awk '/safqmc:/{print $5}' $OUT_FILE >> $VAL_FILE
}