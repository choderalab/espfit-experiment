#!/bin/bash


# Settings
targets="tyk2 cdk2 mcl1 p38"
suffix="spice2-default"
script_path="/home/takabak/data/espfit-experiment/scripts/pl-benchmark"
python ${script_path}/plotall.py --input_prefix "../" --targets "${targets}" --suffix ${suffix}
