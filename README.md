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
- Orientation of the LUT are zero vectors.
- Only supports scanner with one type of module.
- Energy bins are not supported.
- Currently, all event are setted at 1 ms.
- Crystal size is unknown, thus it cannot be used for multiSiddon or Distance Driven projectors.

# TODOs
- Clean up the notation.
- Supports scanner with multiple type of modules.
- Add a check for the number of energy bins.
- Make sure that DetectionBins are not used as CASToR detector elements
- Suport time paquets.
- Enable user to give the scanner a name since many examples do not have scanner name
