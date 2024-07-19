# flow2supera
This repository contains code to translate the HDF5 files output by [ndlar_flow](https://github.com/DUNE/ndlar_flow) to [Supera](https://github.com/DeepLearnPhysics/SuperaAtomic) format for use by the DUNE machine learning reconstruction chain, [lartpc_mlreco3d](https://github.com/DeepLearnPhysics/lartpc_mlreco3d). 

# Prerequisites 

`flow2supera` depends on [edep2supera](https://github.com/DeepLearnPhysics/edep2supera), [SuperaAtomic](https://github.com/DeepLearnPhysics/SuperaAtomic), [larcv](https://github.com/DeepLearnPhysics/larcv2) and [h5flow](https://github.com/peter-madigan/h5flow). Install each of those repositories using the instructions on their respective READMEs and ensure that you can import them in python. Make sure the installation follows this order: `larcv` -> `SuperaAtomic` -> `edep2supera` -> `flow2supera`.

# Installation
Once the prerequisites are met, simply run this command from the top directory:
```
python3 -m pip install .
```

# Usage

The main executable script is located at `bin/run_flow2supera.py` relative to the top directory. The _required_ arguments are the input and output file names and the configuration:
```
python3 bin/run_flow2supera.py -o <output_file> -c 2x2 <input_ndlar_flow_file>
```
Configuration keyword or a file path (full or relative including the file name). Supported configurations: `2x2`, `2x2_data`, `mod1_data`, `2x2_mpvmpr`.
You can also specify the following _optional_ arguments:
- `-n` or `--num_events`: Number of events to process.
- `-s` or `--skip`: Number of first events to skip.
- `-l` or `--log`: Name of a log file to be created.

Upon successful completion, this will produce an output larcv-format file that can be used as input to the machine learning reconstruction. 

# Contributing

Please read the contributing.md file for information on how you can contribute.

# License

Distributed under the MIT License. See LICENSE for more information.
