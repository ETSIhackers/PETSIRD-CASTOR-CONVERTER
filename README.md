# Background
The [Emission Tomography Standardization Initiative (ETSI)](https://etsinitiative.org/)
is working towards establishing a standard for PET Raw Data, called PETSIRD ("PET ETSI Raw Data").
More information is on https://github.com/ETSInitiative/PETSIRD.

# Description
This project provides a converter, in c++, that reads a PETSIRD data file and generates the corresponding CASToR list mode file and geometry LUT.

# Compiling the software
You need to run `yardl generate` in the `PETSIRD/model` directory first.

Compiling the converter:
   ```sh
   cd cpp
   mkdir -p build && cd build`
   cmake -G Ninja -S .. -DHDF5_ROOT=$CONDA_PREFIX
   ninja`
   ```
   If you did not use `conda` to install HDF5, do not add the `-DHDF5_ROOT=$CONDA_PREFIX` part of the `cmake` line.

# How to use
After compilation, run `petsird_castor_converter --help` to see how to use the converter.

# Known limitations
- Multi-TOF kernels or bins are not supported.

# TODOs in the contribution
zxc to fix after hackathon
