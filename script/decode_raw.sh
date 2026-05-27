#!/bin/bash

BIN_NAME=$HOME/spadi/nestdaq_old_type_decoder/install/bin/decode
DATA_PATH=/data/2606_jparc
IN_PATH=$DATA_PATH/raw
OUT_PATH=$DATA_PATH/decode

if [ $# != 1 ]; then
  echo "./shell [run ID]"
  exit 1
fi
run_id=$1

run_id_str=$(printf '%06d' "$run_id")
in_name=$IN_PATH/run$run_id_str.dat
out_name=$OUT_PATH/run${run_id_str}_raw.root

$BIN_NAME $in_name $out_name
