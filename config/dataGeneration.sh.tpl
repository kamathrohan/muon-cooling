#!/bin/bash

if ! echo "$LD_LIBRARY_PATH" | grep -q "/cvmfs/"; then
    source {{ env_setup }}
    echo "Environment Setup!"
else
    echo "Some CVMFS libs are already sourced"
fi

if ! echo "$LD_LIBRARY_PATH" | grep -q "bdsim-install/bin/../lib"; then
    source {{ bdsim_setup }}
    echo "BDSIM Setup!"
else
    echo "BDSIM library path is already in LD_LIBRARY_PATH"
fi

cd $OUTPUT_DIR
bdsim --file=$infile --batch --ngenerate=$ngenerate 
