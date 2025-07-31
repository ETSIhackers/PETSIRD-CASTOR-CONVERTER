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
- CASToR related:
  - Multi-TOF bins are not supported.
  - Energy bins are not supported.
- PETSIRD standard related:
  - Orientation of the LUT are currently defined as zero vectors.
  - Crystal size is unknown, thus it cannot be used for multiSiddon or Distance Driven projectors.
- Others
  - Only supports scanner with one type of module.
  - Currently, all event are setted at 1 ms.

# TODOs
- Soon-ish
  - Clean up the notation.
  - Support time paquets.
  - Crash when detecting multiple types of module.
  - Enable user to give the scanner a name since many examples do not have scanner name.
  - Crash when detecting multiple energy bins.
- Wait for scanner to support these?
  - Supports scanner with multiple type of modules IF they don't differ in detectors size.
  - Make sure that DetectionBins are not used as CASToR detector elements (Needed to support multiple type of modules).
  - Add supports for multi-TOF kernels.
- Future
  - Deduce the size of detectors and propose a voxel size.
