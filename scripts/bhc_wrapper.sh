#!/bin/bash
set -e

#call BHC proper, passing all the arguments
Rscript /scripts/bhc.R ${@:1}

#post-hoc processing to make BiNGO and MEME outputs
mkdir functional_analysis_inputs
python3 /scripts/bingomeme.py