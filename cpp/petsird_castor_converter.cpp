#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include "petsird_helpers.h"
#include "petsird_helpers/geometry.h"
#include "generated/hdf5/protocols.h"
#include "generated/types.h"
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
//#include <xtensor/xview.hpp>

#include "generated/binary/protocols.h"
using petsird::binary::PETSIRDReader;
using namespace std;

// // ===============================================================
// // ===============================================================
// // Function: get_detection_efficiency()
// // Arguments: ScannerInformation, uint32_t, uint32_t
// //                              , uint32_t, uint32_t
// // Return: float
// // Description: This function overloads the one from helpers
// //              by simply using the two uint32_t detector indices
// //              supplied in arguments to retrieve the global
// //              detection efficiency of the LOR and return it.
// //              It also uses the energy indices supplied.
// // ===============================================================
// // ===============================================================
// float get_detection_efficiency( const petsird::ScannerInformation& scanner, uint32_t detector_id1, uint32_t detector_id2,
//                                                                             uint32_t energy_id1,   uint32_t energy_id2 )
// {
//   // The global efficiency of the detector pair
//   float eff = 1.0F;
//   // Get the detector element efficiencies from the scanner
//   const auto& det_el_efficiencies = scanner.detection_efficiencies.det_el_efficiencies;
//   if (det_el_efficiencies)
//   {
//     // Multiply global efficiency of the pair by the individual efficiency of each detector
//     eff *= ( (*det_el_efficiencies)(detector_id1, energy_id1)
//            * (*det_el_efficiencies)(detector_id2, energy_id2) );
//   }
//   // Get the module pair efficiencies from the scanner
//   const auto& module_pair_efficiencies_vector = scanner.detection_efficiencies.module_pair_efficiencies_vector;
//   if (module_pair_efficiencies_vector)
//   {
//     // Get the SGID LUT from the scanner
//     const auto& module_pair_SGID_LUT = scanner.detection_efficiencies.module_pair_sgidlut;
//     assert(module_pair_SGID_LUT);
//     // Get module and element ids for each detector in the pair
//     vector<uint32_t> detector_ids;
//     detector_ids.push_back(detector_id1);
//     detector_ids.push_back(detector_id2);
//     const auto mod_and_els = petsird_helpers::get_module_and_element(scanner.scanner_geometry, detector_ids);
//     assert(scanner.scanner_geometry.replicated_modules.size() == 1);
//     // Get the SGID of this module pair
//     const int SGID = (*module_pair_SGID_LUT)(mod_and_els[0].module, mod_and_els[1].module);
//     // When SGID is negative, this means that the pair of detectors cannot be in coincidence
//     if (SGID<0) return 0.;
//     // Otherwise, retrieve module pair efficiencies for this SGID
//     const auto& module_pair_efficiencies = (*module_pair_efficiencies_vector)[SGID];
//     assert(module_pair_efficiencies.sgid == static_cast<unsigned>(SGID));
//     // Multiply the global efficiency by the efficiency of the module pair for this pair of detector elements
//     eff *= module_pair_efficiencies.values(mod_and_els[0].el, energy_id1, mod_and_els[1].el, energy_id2);
//   }
//   return eff;
// }
// // ===============================================================
// // ===============================================================
// // Function: transform_Coordinate()
// // Arguments: Coordinate, RigidTransformation, RigidTransformation
// // Return: Coordinate
// // Description: The two supplied rigid transformations are applied
// //              sequentially on the supplied coordinate and the
// //              resulting transformed coordinate is returned.
// // ===============================================================
// // ===============================================================
// petsird::Coordinate transform_Coordinate(const petsird::Coordinate& coord, const petsird::RigidTransformation& volume_transform, const petsird::RigidTransformation& module_transform)
// {
//     // Buffer coordinates
//     petsird::Coordinate transformed_coord1;
//     petsird::Coordinate transformed_coord2;
//     // Apply volume rotation
//     for (uint32_t i = 0; i < 3; i++)
//     {
//         // Reset coord
//         transformed_coord1.c[i] = 0.;
//         // Apply rotation
//         for (uint32_t j = 0; j < 3; j++)
//             transformed_coord1.c[i] += volume_transform.matrix(i, j) * coord.c[j];
//     }
//     // Apply volume translation
//     for (uint32_t i = 0; i < 3; i++)
//     {
//         transformed_coord1.c[i] += volume_transform.matrix(i, 3);
//     }
//     // Apply module rotation
//     for (uint32_t i = 0; i < 3; i++)
//     {
//         // Reset coord
//         transformed_coord2.c[i] = 0.;
//         // Apply rotation
//         for (uint32_t j = 0; j < 3; j++)
//             transformed_coord2.c[i] += module_transform.matrix(i, j) * transformed_coord1.c[j];
//     }
//     // Apply module translation
//     for (uint32_t i = 0; i < 3; i++)
//     {
//         transformed_coord2.c[i] += module_transform.matrix(i, 3);
//     }
//     // Return transformed coordinate
//     return transformed_coord2;
// }
// // ===============================================================
// // ===============================================================
// // Function: transform_BoxShape()
// // Arguments: BoxShape, RigidTransformation, RigidTransformation
// // Return: BoxShape
// // Description: The two supplied rigid transformations are applied
// //              on each of the 8 coordinates of the supplied
// //              box shape using the transfor_Coordinate function
// //              and the resulting transformed box shape is
// //              returned.
// // ===============================================================
// // ===============================================================
// petsird::BoxShape transform_BoxShape(const petsird::BoxShape& box_shape, const petsird::RigidTransformation& volume_transform, const petsird::RigidTransformation& module_transform)
// {
//     // The transformed box shape
//     petsird::BoxShape transformed_box;
//     // Loop on the 8 box corners
//     for (int i=0; i<8; i++)
//     {
//         // Apply the transformation on each corner
//         transformed_box.corners[i] = transform_Coordinate(box_shape.corners[i], volume_transform, module_transform);
//     }
//     // Return the transformed box shape
//     return transformed_box;
// }
// // Temporary function to apply a manual shift on a box shape
// petsird::BoxShape translate_BoxShape(const petsird::BoxShape& box_shape)
// {
//     petsird::BoxShape translated_box;
//     for (int i=0; i<8; i++)
//     {
//       petsird::Coordinate transformed_coord1;
//       transformed_coord1.c[0] = box_shape.corners[i].c[0];
//       transformed_coord1.c[1] = box_shape.corners[i].c[1] + 1.6;
//       transformed_coord1.c[2] = box_shape.corners[i].c[2] + 1.6;
//       //cout << "coord orig: " << box_shape.corners[i].c[1] << " | translated coord: " << transformed_coord1.c[1] << endl;
//       translated_box.corners[i] = transformed_coord1;
//     }
//     //cout << "translated" << endl << flush;
//     return translated_box;
// }
// // ===============================================================
// ===============================================================
// Function: compute_centroid_BoxShape()
// Arguments: BoxShape
// Return: Coordinate
// Description: The centroid of the 8 corners of the supplied box
//              shape is computed and returned as a coordinate.
// ===============================================================
// ===============================================================
petsird::Coordinate compute_centroid_BoxShape(const petsird::BoxShape& box_shape)
{
//  cout << "=============== COMPUTE CENTROID" << endl;
    // The centroid coordinate
    petsird::Coordinate centroid_coord;
    // Loop on the coordinates
    for (int i = 0; i < 3; i++)
    {
        // Reset coordinate
        centroid_coord.c[i] = 0.;
        // Loop on the 8 box corners
        for (int c = 0; c < 8; c++)
        {
            // Add corner contribution
            centroid_coord.c[i] += (box_shape.corners[c]).c[i];
//            if (i==0) cout << "coord[" << c << "]  = " << (box_shape.corners[c]).c[i] << endl;
        }
        // Average
        centroid_coord.c[i] /= 8.;
//        if (i==0) cout << "centroid: " << centroid_coord.c[i] << endl;
    }
    // Return centroid
    return centroid_coord;
}
// ===============================================================
// ===============================================================
// Function: show_help()
// Description: Prints out help about usage and option on screen.
// ===============================================================
// ===============================================================
void show_help()
{
  cout << endl;
  cout << "Usage:  petsird_castor_converter  -in input_file.petsird  [-out output_base_name]  [options]" << endl;
  cout << endl;
  cout << "[Description]:" << endl;
  cout << "  This program can be used to convert a PETSIRD file into files suitable for CASToR." << endl;
  cout << "  From the PETSIRD file, the program can output a CASToR scanner LUT description that" << endl;
  cout << "  follows the scanner description embedded into the PETSIRD file, a CASToR normalization" << endl;
  cout << "  datafile describing the set of all detector pairs in coincidence with their respective" << endl;
  cout << "  efficiencies, and a CASToR listmode datafile containing all prompt events from the" << endl;
  cout << "  PETSIRD file. By default, the program only goes through the scanner description. Flags" << endl;
  cout << "  must be used to enable actions and writing of the different output. Otherwise, verbose" << endl;
  cout << "  levels can be upped to print detailed information of the content of the PETSIRD file." << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -in file.petsird     : Provide the input petsird file to be converter" << endl;
  cout << "  -out output          : Provide an output base name used as a prefix for all output files" << endl;
  cout << "                         (needed with option -write)." << endl;
  cout << "  -norm                : Go into normalization/efficiency description." << endl;
  cout << "  -list                : Go into listmode data." << endl;
  cout << "  -write               : Actually write the CASToR output associated to each action." << endl;
  cout << "  -vb level            : Give an integer verbose level; 0 nothing, 1 global, 2 more details," << endl;
  cout << "                         3 printing in loops, 4 even more printing, 5 become interactive (default: 1)." << endl;
  cout << "  -h, --help           : Prints this help." << endl;
  cout << endl;
}
// ======================================================================================================================
// ======================================================================================================================
// ======================================================================================================================
// M A I N   P R O G R A M
// ======================================================================================================================
// ======================================================================================================================
// ======================================================================================================================
int main(int argc, char** argv)
{
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------

  // Show help if no argument
  if (argc==1)
  {
    show_help();
    return 0;
  }

  // ----------------------------------
  // Lists of parameters
  // ----------------------------------

  // Verbose level (default is 1)
  int verbose = 1;
  // Input petsird file
  string petsird_file = "";
  // Output base name
  string output_base_name = "";
  // Flag to go into normalization/efficiency
  bool flag_norm = false;
  // Flag to go into listmode data
  bool flag_list = false;
  // Flag to say that we write some output
  bool flag_write = false;

  // ----------------------------------
  // Read parameters list
  // ----------------------------------

  for (int i=1; i<argc; i++)
  {
    // Get option
    string option = (string)argv[i];
    // Input PETSIRD file
    if (option=="-in")
    {
      if (i==argc-1)
      {
        cerr << "***** Argument missing after option '" << option << "' !" << endl;
        exit(1);
      }
      petsird_file = (string)argv[i+1];
      i++;
    }
    // Output base name
    else if (option=="-out")
    {
      if (i==argc-1)
      {
        cerr << "***** Argument missing after option '" << option << "' !" << endl;
        exit(1);
      }
      output_base_name = (string)argv[i+1];
      i++;
    }
    // Flags to enable actions
    else if (option=="-norm") flag_norm = true;
    else if (option=="-list") flag_list = true;
    else if (option=="-write") flag_write = true;
    // Verbose level
    else if (option=="-vb")
    {
      if (i==argc-1)
      {
        cerr << "***** Argument missing after option '" << option << "' !" << endl;
        exit(1);
      }
      verbose = atoi(argv[i+1]);
      i++;
    }
    // Show help
    else if (option=="-h" || option=="--help" || option=="--h" || option=="-help")
    {
      show_help();
      return 0;
    }
    // Unknown option
    else
    {
      cerr << "***** Unknown option '" << option << "' !" << endl;
      cerr << "      Please see help using -h option." << endl;
      exit(1);
    }
  }

  // ----------------------------------
  // Checks about options
  // ----------------------------------

  // Required input PETSIRD file
  if (petsird_file=="")
  {
    cerr << "***** An input PETSIRD file is required !" << endl;
    exit(1);
  }
  // Required output base name
  if (flag_write && output_base_name=="")
  {
    cerr << "***** An output base name is required !" << endl;
    exit(1);
  }
  // Verbose level
  if (verbose<0)
  {
    cerr << "***** Positive verbose level required !" << endl;
    exit(1);
  }

  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------

  // Verbose
  if (verbose>1) cout << "====================================================" << endl;
  if (verbose>0) cout << "==  Initializing PETSIRD file" << endl;
  if (verbose>1) cout << "====================================================" << endl;
  if (verbose>1) cout << "--> Input PETSIRD file: " << petsird_file << endl;
  // Open the file
  PETSIRDReader reader(petsird_file.c_str());
  // Read header
  petsird::Header header;
  reader.ReadHeader(header);

  // Get scanner name
  string scanner_name = header.scanner.model_name;
  // If scanner name is empty
  if (scanner_name=="")
  {
    // If not writing things, then use a generic one
    if (!flag_write) scanner_name = "petsird_scanner";
    // If writing things, then use the output base name
    else scanner_name = output_base_name;
  }
  if (verbose>1) cout << "--> Scanner model: " << scanner_name << endl;

  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------

  // Action
  if (verbose>1) cout << "====================================================" << endl;
  if (verbose>0) cout << "==  Going through scanner description" << endl;
  if (verbose>1) cout << "====================================================" << endl;

  // CASToR LUT scanner file
  uint32_t nb_data_written_in_lut = 0;
  FILE* flut = NULL;
  if (flag_write)
  {
    string scanner_lut_file = scanner_name + ".lut";
    flut = fopen(scanner_lut_file.c_str(), "wb");
    if (flut==NULL)
    {
      cerr << "***** Failed to create output CASToR LUT file '" << scanner_lut_file << "' !" << endl;
      exit(1);
    }
  }

  /*
  // Need a zero coordinate
  // zxc not used?
  petsird::Coordinate zero_coord;
  zero_coord.c[0] = 0.;
  zero_coord.c[1] = 0.;
  zero_coord.c[2] = 0.;
  */

  // Declare the maps to get from petsird 4 IDs to CASToR flatenned ID
  // zxc Now it is three? Type Module, module instances, detector instance?
  map<tuple<petsird::TypeOfModule, petsird::DetectionBin>, uint32_t> map_petsird2castor_id;
  map<uint32_t, tuple<petsird::TypeOfModule, petsird::DetectionBin>> map_castor2petsird_id;

  // The table to get castor ids from petsird ids
  int nb_detectors_in_scanner = 0;
  for (petsird::TypeOfModule tm=0; tm<header.scanner.scanner_geometry.NumberOfReplicatedModules(); tm++)
  {
    nb_detectors_in_scanner += petsird_helpers::get_num_det_els(header.scanner, tm);
  }
  if (verbose>1) cout << "--> Number of detector in scanner : " << nb_detectors_in_scanner << endl;
  //OBS uint32_t**** table_petsird2castor_id = (uint32_t****)malloc(nb_detectors_in_scanner * sizeof(uint32_t***));

  // We will search for the minimum and maximum positions in each axis independently
  // We initialize the values with the centroid position of the first detecting element
  const petsird::TypeOfModule type_of_module_zero{ 0 };
  const petsird::DetectionBin detection_bin_zero{ 0 };
  auto coord_init = compute_centroid_BoxShape(petsird_helpers::geometry::get_detecting_box(header.scanner, type_of_module_zero, detection_bin_zero));
  float min_x = coord_init.c[0];
  float min_y = coord_init.c[1];
  float min_z = coord_init.c[2];
  float max_x = coord_init.c[0];
  float max_y = coord_init.c[1];
  float max_z = coord_init.c[2];
  // We will compute the center of mass of the scanner
  double center_x = 0.;
  double center_y = 0.;
  double center_z = 0.;

  // Loop over the lists of replicates module (top level in the scanner geomety)
  if (verbose>1) cout << "--> Processing scanner elements" << endl;
  uint32_t flat_index = 0;
  petsird::ExpandedDetectionBin expanded_detection_bin;
  // Unsure if init is needed
  expanded_detection_bin.energy_index = 0;
  for (uint32_t typeMod=0; typeMod<header.scanner.scanner_geometry.NumberOfReplicatedModules(); typeMod++)
  {
    auto& rep_module = header.scanner.scanner_geometry.replicated_modules[typeMod];
    auto& det_els = rep_module.object.detecting_elements;
    for (uint32_t repMod=0; repMod<rep_module.transforms.size(); repMod++)
    {
      expanded_detection_bin.module_index = repMod;
      for (uint32_t detEl=0; detEl<det_els.transforms.size(); detEl++)
      {
        expanded_detection_bin.element_index = detEl;
        // zxc typeMod Should be object? Seems ok in ana
        auto box = petsird_helpers::geometry::get_detecting_box(header.scanner, typeMod, expanded_detection_bin);
        // Compute centroid
        auto centroid = compute_centroid_BoxShape(box);
        // Verbose
        if (verbose>2)
        {
          cout << "          | Detector id: " << flat_index << " | PETSIRD detector ids [ " << typeMod << " ; " << repMod << " ; " << detEl << " ]" << endl;
          cout << "          | Centroid [" << centroid.c[0] << "; " << centroid.c[1] << "; " << centroid.c[2] << "] mm" << endl;
        }
        if (verbose>4)
        {
          cout << "          | Press enter to continue" << endl;
          getchar();
        }
        // Write centroid into the CASToR LUT file
        if (flag_write)
        {
          // This is the centroid of the detecting element
          nb_data_written_in_lut += fwrite(&(centroid.c[0]), sizeof(float), 1, flut);
          nb_data_written_in_lut += fwrite(&(centroid.c[1]), sizeof(float), 1, flut);
          nb_data_written_in_lut += fwrite(&(centroid.c[2]), sizeof(float), 1, flut);
          // The orientation vector is set to zero for now
          float zero = 0.;
          nb_data_written_in_lut += fwrite(&zero, sizeof(float), 1, flut);
          nb_data_written_in_lut += fwrite(&zero, sizeof(float), 1, flut);
          nb_data_written_in_lut += fwrite(&zero, sizeof(float), 1, flut);
        }
        // Add the indices to the map
        petsird::DetectionBin detection_bin = petsird_helpers::make_detection_bin(header.scanner, typeMod, expanded_detection_bin);
        tuple<petsird::TypeOfModule, petsird::DetectionBin> petsird_ids = tuple(typeMod,detection_bin);
        map_petsird2castor_id.insert({ petsird_ids, flat_index });
        map_castor2petsird_id.insert({ flat_index, petsird_ids });
        //OBS table_petsird2castor_id[ml][mtr][dl][dtr] = flat_index;
        // Increment the flat index
        flat_index++;
        // Update min/max x,y,z
        if (centroid.c[0]<min_x) min_x = centroid.c[0];
        if (centroid.c[1]<min_y) min_y = centroid.c[1];
        if (centroid.c[2]<min_z) min_z = centroid.c[2];
        if (centroid.c[0]>max_x) max_x = centroid.c[0];
        if (centroid.c[1]>max_y) max_y = centroid.c[1];
        if (centroid.c[2]>max_z) max_z = centroid.c[2];
        // Update center of mass
        center_x += ((double)(centroid.c[0]));
        center_y += ((double)(centroid.c[1]));
        center_z += ((double)(centroid.c[2]));
      }
    }
  }

  // Finish center of mass computation
  center_x /= ((double)flat_index);
  center_y /= ((double)flat_index);
  center_z /= ((double)flat_index);

  // Verbose
  if (verbose>1)
  {
    cout << "--> Number of detecting elements: " << flat_index << endl;
    cout << "--> Min/Max of detector elements centroid coordinates" << endl;
    cout << "  | Along X axis [ " << min_x << " ; " << max_x << " ]" << endl;
    cout << "  | Along Y axis [ " << min_y << " ; " << max_y << " ]" << endl;
    cout << "  | Along Z axis [ " << min_z << " ; " << max_z << " ]" << endl;
    cout << "--> Scanner center of mass: [ " << center_x << " ; " << center_y << " ; " << center_z << " ]" << endl;
  }

  // Close scanner LUT file
  if (flag_write)
  {
    fclose(flut);
    // Check for accurate writing
    if (nb_data_written_in_lut != flat_index * 6)
    {
        cerr << "***** Failed to write all data in scanner LUT file !" << endl;
        return 1;
    }
    // Verbose
    if (verbose>1) cout << "--> CASToR scanner LUT file '" << scanner_name << ".lut' written" << endl;
  }
  // Write scanner header
  if (flag_write)
  {
    string scanner_hscan_file = scanner_name + ".hscan";
    ofstream hscan(scanner_hscan_file.c_str());
    if (!hscan)
    {
      cerr << "***** Failed to open scanner header file '" << scanner_hscan_file << "' for writing !" << endl;
      exit(1);
    }
    hscan << "modality: PET" << endl;
    hscan << "number of elements: " << flat_index << endl;
    hscan << "scanner radius: " << max(max_y-min_y,max_x-min_x) << endl;
    hscan << "voxels number transaxial: 100" << endl;
    hscan << "voxels number axial: 100" << endl;
    hscan << "field of view transaxial: " << max(max_y-min_y,max_x-min_x) << endl;
    hscan << "field of view axial: " << max_z-min_z << endl;
    hscan << "mean depth of interaction: -1" << endl;
    hscan << "######################################################" << endl;
    hscan << "# The following lines are irrelevant because cannot be" << endl;
    hscan << "# known from the petsird file but required by castor" << endl;
    hscan << "######################################################" << endl;
    hscan << "crystals size depth: 0." << endl;
    hscan << "number of rings in scanner: 1" << endl;
    hscan << "number of layers: 1" << endl;
    hscan << "number of crystals in layer: " << flat_index << endl;
    hscan.close();

    // Verbose
    if (verbose>1) cout << "--> CASToR scanner header file '" << scanner_hscan_file << "' written" << endl;
  }

//   // -------------------------------------------------------------------------------
//   // -------------------------------------------------------------------------------
//   // -------------------------------------------------------------------------------
//   // -------------------------------------------------------------------------------
//   // -------------------------------------------------------------------------------
//   // -------------------------------------------------------------------------------

//   // Flag for going into normalization components
//   if (flag_norm)
//   {
//     // Verbose
//     if (verbose>1) cout << "====================================================" << endl;
//     if (verbose>0) cout << "==  Going through normalization components" << endl;
//     if (verbose>1) cout << "====================================================" << endl;

//     // Open castor normalization file
//     string castor_norm_file = output_base_name + "_norm.cdf";
//     FILE* fnorm = NULL;
//     if (flag_write)
//     {
//       fnorm = fopen(castor_norm_file.c_str(), "wb");
//       if (fnorm==NULL)
//       {
//         cerr << "***** Failed to create output CASToR normalization file '" << castor_norm_file << "' !" << endl;
//         exit(1);
//       }
//     }
//     uint64_t nb_data_written_in_norm = 0;
//     uint64_t nb_valid_lors = 0;
//     uint64_t size_of_a_norm_lor = 3;

//     // Start a double loop over the detector elements with duplicates
//     if (verbose>1) cout << "--> Loop over all possible detector elements pairs" << endl;
//     for (uint32_t id1=0; id1<flat_index; id1++)
//     {
//       for (uint32_t id2=id1+1; id2<flat_index; id2++)
//       {
//         // Get the detection efficiency for this detector pair (TODO: manage energies)
//         float eff = get_detection_efficiency( header.scanner, id1, id2, 0, 0 );
//         // Verbose
//         if (verbose>2) cout << "  --> Pair [ " << id1 << " ; " << id2 << " ] efficiency is " << eff << endl << flush;
// 	if (verbose>4)
// 	{
//           cout << "    | Press enter to continue" << endl;
//           getchar();
// 	}
//         // Consider the LOR valid only if efficiency is strictly positive
//         if (eff>0.)
//         {
//           // Compute normalization factor as the inverse of efficiency
//           float norm = 1. / eff;
//           // Write data
// 	  if (flag_write)
// 	  {
//             nb_data_written_in_norm += fwrite(&norm,sizeof(float),1,fnorm);
//             nb_data_written_in_norm += fwrite(&id1,sizeof(uint32_t),1,fnorm);
//             nb_data_written_in_norm += fwrite(&id2,sizeof(uint32_t),1,fnorm);
// 	  }
//           // Increment the number of valid lors
//           nb_valid_lors++;
//         }
//       }
//     }

//     // Verbose
//     if (verbose>1) cout << "--> Number of valid LORs: " << nb_valid_lors << endl;

//     // Close CASToR norm file
//     if (flag_write)
//     {
//       fclose(fnorm);
//       if (nb_data_written_in_norm!=nb_valid_lors*size_of_a_norm_lor)
//       {
//         cerr << "***** Failed to write all data into CASToR normalization file !" << endl;
//         exit(1);
//       }
//       if (verbose>1) cout << "--> CASToR normalization file '" << castor_norm_file << "' written" << endl;
//     }

//     // Write castor header file
//     if (flag_write)
//     {
//       string castor_hnorm_file = output_base_name + "_norm.cdh";
//       ofstream cdhn(castor_hnorm_file.c_str());
//       if (!cdhn)
//       {
//         cerr << "***** Failed to create castor normalization header file '" << castor_hnorm_file << "' for writing !" << endl;
//         exit(1);
//       }
//       cdhn << "Scanner name: " << scanner_name << endl;
//       cdhn << "Data filename: " << castor_norm_file << endl;
//       cdhn << "Number of events: " << nb_valid_lors << endl;
//       cdhn << "Data mode: normalization" << endl;
//       cdhn << "Data type: PET" << endl;
//       cdhn << "Start time (s): 0" << endl;
//       cdhn << "Duration (s): 1" << endl;
//       cdhn << "Normalization correction flag: 1" << endl;
//       cdhn.close();
//       // Verbose
//       if (verbose>1) cout << "--> CASToR normalization header file '" << castor_hnorm_file << "' written" << endl;
//     }

//   // End if going through normalization components
//   }

  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------

  // Flag for going through listmode data
  if (flag_list)
  {
    // Verbose
    if (verbose>1) cout << "====================================================" << endl;
    if (verbose>0) cout << "==  Processing PETSIRD listmode data" << endl;
    if (verbose>1) cout << "====================================================" << endl;

    // ---------------------------
    // TOF management
    // ---------------------------
    bool tof_enable = false;
    double speed_of_light_in_mm_per_ps = 0.299792458;
    double tof_measurement_range_in_ps = -1.;
    double tof_resolution_in_ps = -1.;
    bool tof_quantization = true;
    double tof_quantization_bin_size_in_ps = -1.;
    // Search for the higher number of TOF bins through pairs of module types
    int header_scanner_NumberOfTOFBins = 1;
    bool varying_number_of_tof_bins = false;
    for (int rm1=0; rm1<header.scanner.scanner_geometry.NumberOfReplicatedModules(); rm1++)
    {
      for (int rm2=rm1; rm2<header.scanner.scanner_geometry.NumberOfReplicatedModules(); rm2++)
      {
        int current_NumberOfTOFBins = header.scanner.tof_bin_edges[rm1][rm2].NumberOfBins();
        if (current_NumberOfTOFBins > header_scanner_NumberOfTOFBins)
        {
          cout << "rm1: " << rm1 << " | rm2: " << rm2 << " | number of TOF bins: " << header.scanner.tof_bin_edges[rm1][rm2].NumberOfBins() << endl;
          if (header_scanner_NumberOfTOFBins!=1) varying_number_of_tof_bins = true;
          header_scanner_NumberOfTOFBins = current_NumberOfTOFBins;
        }
      }
    }

    // If more than one TOF bin, then we assume TOF is enable
    if (header_scanner_NumberOfTOFBins!=1)
    {
      // IMPORTANT NOTE: CASToR works in the following way wrt TOF. For histogram/sinogram,
      // the TOF bin size is unique. For listmode data, the measurement is supposed to be a
      // value which can or cannot by quantized. If quantized, then, the quantization bin
      // size can be supplied in the header of the datafile which can be taken into account
      // in the projector. For the conversion of PETSIRD to CASToR data, we will use this
      // quantization bin size as the PETSIRD bin size. Then, each event will be affected
      // by a value at the middle of the associated TOF bin. The quantization bin size is
      // supposed to be fixed in CASToR, so we check for the standard deviation of the
      // PETSIRD TOF bin sizes to decide if the hypothesis holds or not, based on a
      // threshold value that can be adjusted below. This threshold should not be zero to
      // avoid numerical rounding problems. In case the PETSIRD data use varying TOF bin
      // size, a workaround could be used by not providing the quantization bin size in
      // CASToR. That way, CASToR will assume that TOF measurements provided for each event
      // are continuous, so the projector will simply use a Gaussian function centered on
      // the measurement to model the TOF function. This is what is implemented below but
      // we still throw a warning to warn the user.
      //
      // TO DO: for the moment the following code only works for module types that have the
      //        same TOF resolution and number of bins. So we crash if not the case.
      if (varying_number_of_tof_bins)
      {
        cerr << "***** Multiple TOF number of bins not supported YET !" << endl;
        return 1;
      }
      float tof_resolution_default = header.scanner.tof_resolution[0][0];
      if (header.scanner.scanner_geometry.NumberOfReplicatedModules()>1)
      {
        for (int rm1=1; rm1<header.scanner.scanner_geometry.NumberOfReplicatedModules(); rm1++)
        {
          for (int rm2=rm1; rm2<header.scanner.scanner_geometry.NumberOfReplicatedModules(); rm2++)
          {
            if (header.scanner.tof_resolution[rm1][rm2]!=tof_resolution_default)
            {
              cerr << "***** Varying TOF resolution not supported, YET !" << endl;
              return 1;
            }
          }
        }
      }
      tof_enable = true;
      // TOF resolution in ps (needed by castor)
      tof_resolution_in_ps = tof_resolution_default * 2. / speed_of_light_in_mm_per_ps;
      // Verbose
      if (verbose>1)
      {
        cout << "--> TOF FWHM resolution: " << tof_resolution_default << " mm (" << tof_resolution_in_ps << " ps)" << endl;
        cout << "--> Number of TOF bins: " << header_scanner_NumberOfTOFBins << endl;
        cout << "--> PETSIRD TOF bin edges: " << header.scanner.tof_bin_edges[0][0].edges << endl;
      }
      // Compute mean TOF bin size in mm from petsird header values
      double mean_tof_bin_size_in_mm = 0.;
      for (uint32_t tb=0; tb<header_scanner_NumberOfTOFBins; tb++)
        mean_tof_bin_size_in_mm += header.scanner.tof_bin_edges[0][0].edges[tb+1] - header.scanner.tof_bin_edges[0][0].edges[tb];
      mean_tof_bin_size_in_mm /= ((double)(header_scanner_NumberOfTOFBins));
      // Compute standard deviation of the TOF bin size
      double stdv_tof_bin_size_in_mm = 0.;
      for (uint32_t tb=0; tb<header_scanner_NumberOfTOFBins; tb++)
        stdv_tof_bin_size_in_mm += pow((mean_tof_bin_size_in_mm - header.scanner.tof_bin_edges[0][0].edges[tb+1] + header.scanner.tof_bin_edges[0][0].edges[tb]) , 2.);
      stdv_tof_bin_size_in_mm /= ((double)(header_scanner_NumberOfTOFBins));
      // Check hypothesis that the TOF bin size is fixed
      double difference_in_tof_bin_size_tolerance_in_percent = 0.1;
      if (100. * stdv_tof_bin_size_in_mm / mean_tof_bin_size_in_mm > difference_in_tof_bin_size_tolerance_in_percent)
      {
        // TOF bin size not fixed: see note above
        cerr << "!!!!! TOF bin size of PETSIRD data has been detected to be variable whereas CASToR currently only works with fixed bin size." << endl;
        cerr << "!!!!! The measurement of the TOF for each event will thus be assumed to be continuous and the value provided will be the one" << endl;
        cerr << "!!!!! at the center of the TOF bin. For more information, see the CASToR TOF related documentation or the technical note by" << endl;
        cerr << "!!!!! Filipovic et al, Phys Med Biol, 2019: 'Time-of-flight (TOF) implementation for PET reconstruction in practice'." << endl;
        // Measurements assumed to be continuous as opposed to quantized
        tof_quantization = false;
      }
      else
      {
        // Affect the CASToR quantization bin size as the mean TOF bin size from PETSIRD data
        tof_quantization_bin_size_in_ps = mean_tof_bin_size_in_mm * 2. / speed_of_light_in_mm_per_ps;
        // Verbose
        if (verbose>1) cout << "--> CASToR quantization bin size: " << tof_quantization_bin_size_in_ps << " ps" << endl;
      }
      // Compute the range of TOF measurements for CASToR
      tof_measurement_range_in_ps = ( header.scanner.tof_bin_edges[0][0].edges[header_scanner_NumberOfTOFBins] - header.scanner.tof_bin_edges[0][0].edges[0] )
                                  * 2. / speed_of_light_in_mm_per_ps;
      // Verbose
      if (verbose>1) cout << "--> CASToR TOF measurement range: " << tof_measurement_range_in_ps << " ps" << endl;
    }

    // Process events by time blocks
    float energy_1 = 0, energy_2 = 0;
    std::size_t num_prompts = 0;
    float last_time = 0.F;
    petsird::TimeBlock time_block;
    uint32_t time_block_index = 0;

    // Open castor data file
    string castor_data_file = output_base_name + ".cdf";
    FILE* fcastor = NULL;
    if (flag_write)
    {
      fcastor = fopen(castor_data_file.c_str(), "wb");
      if (fcastor==NULL)
      {
        cerr << "***** Failed to create output CASToR data file '" << castor_data_file << "' !" << endl;
        exit(1);
      }
    }
    int nb_data_written = 0;
    int nb_events = 0;

    // Loop on time blocks
    if (verbose>1) cout << "--> Processing time blocks" << endl;
    while (reader.ReadTimeBlocks(time_block))
    {
      // Check if it is an actual event time block
      if (std::holds_alternative<petsird::EventTimeBlock>(time_block))
      {
        auto& event_time_block = std::get<petsird::EventTimeBlock>(time_block);
        last_time = event_time_block.time_interval.stop;
        if (verbose>2) cout << "  --> Time block index " << time_block_index << " with time " << last_time << endl << flush;
        const petsird::TypeOfModulePair type_of_module_pair{ 0, 0 }; 
        num_prompts += event_time_block.prompt_events[type_of_module_pair[0]][type_of_module_pair[1]].size(); 
        if (verbose>=5) 
            std::cout << "=====================  Prompt events in time block from " << last_time << " ==============\n"; 
        if (verbose>2) cout << "    | Number of prompts: " << num_prompts << endl << flush;
        if (verbose>5)
        {
          cout << "    | Press enter to continue" << endl;
          getchar();
        }
        // Note: just doing one module-type ATM 
        for (auto& event : event_time_block.prompt_events[type_of_module_pair[0]][type_of_module_pair[1]]) 
        { 
           const auto expanded_detection_bin0 
                = petsird_helpers::expand_detection_bin(header.scanner, type_of_module_pair[0], event.detection_bins[0]); 
           const auto expanded_detection_bin1 
                = petsird_helpers::expand_detection_bin(header.scanner, type_of_module_pair[1], event.detection_bins[1]); 
//TODO           energy_1 += energy_mid_points[expanded_detection_bin0.energy_index]; 
//TODO           energy_2 += energy_mid_points[expanded_detection_bin1.energy_index]; 

//           //std::cout << "CoincidenceEvent(detectorIds=[" << event.detection_bins[0] << ", " << event.detection_bins[1] << "], tofIdx=" << event.tof_idx << "])\n";
//           const auto module_and_elems
//               = petsird_helpers::get_module_and_element(header.scanner.scanner_geometry, event.detector_ids);
//           //std::cout << "    " << "[ModuleAndElement(module=" << module_and_elems[0].module << ", " << "el=" << module_and_elems[0].el << ")," << " ModuleAndElement(module = " << module_and_elems[1].module << ", " << "el=" << module_and_elems[1].el << ")]\n";
//           //std::cout << "    efficiency:" << petsird_helpers::get_detection_efficiency(header.scanner, event) << "\n";

           // Make petsird tuple indices
           tuple<petsird::TypeOfModule, petsird::DetectionBin> petsird_id1 = tuple(0, event.detection_bins[0]);
           tuple<petsird::TypeOfModule, petsird::DetectionBin> petsird_id2 = tuple(0, event.detection_bins[1]);
           auto search1 = map_petsird2castor_id.find(petsird_id1);
           auto search2 = map_petsird2castor_id.find(petsird_id2);
           uint32_t castor_id1 = search1->second;
           uint32_t castor_id2 = search2->second;

//OBS           castor_id1 = table_petsird2castor_id[0][module_and_elems[0].module][0][module_and_elems[0].el];
//OBS           castor_id2 = table_petsird2castor_id[0][module_and_elems[1].module][0][module_and_elems[1].el];
// //          cout << "  castor_id1: " << castor_id1 << " | petsird_id1: " << event.detector_ids[0] << endl;
// //          cout << "  castor_id2: " << castor_id2 << " | petsird_id2: " << event.detector_ids[1] << endl;
// /*
//                 // TOF management
//                 if (tof_enable)
//                 {
//                   uint32_t petsird_tof_bin_id = event.tof_idx;
//                   cout << "tof index: " << petsird_tof_bin_id << endl;
//                   cout << "tof bin edges: [ " << header.scanner.tof_bin_edges[petsird_tof_bin_id] << " ; " << header.scanner.tof_bin_edges[petsird_tof_bin_id+1] << " ]" << endl;
//                   float petsird_tof_middle_point_in_mm = (header.scanner.tof_bin_edges[petsird_tof_bin_id+1] + header.scanner.tof_bin_edges[petsird_tof_bin_id]) / 2.;
//                   float castor_tof_measurement_in_ps = petsird_tof_middle_point_in_mm * 2. / speed_of_light_in_mm_per_ps;
//                   cout << "tof middle point: " << petsird_tof_middle_point_in_mm << " mm" << endl;
//                   cout << "castor tof measurement: " << castor_tof_measurement_in_ps << " ps " << endl;
//                   getchar();
//                 }
//                 */
          // Write a castor event
          if (flag_write)
          {
            uint32_t time_zero = 1;
            // TODO look if time is correct in PETSIRD file and transfer it into castor datafile
            nb_data_written += fwrite(&time_zero, sizeof(uint32_t), 1, fcastor);
            nb_data_written += fwrite(&castor_id1, sizeof(uint32_t), 1, fcastor);
            nb_data_written += fwrite(&castor_id2, sizeof(uint32_t), 1, fcastor);
          }
          // Increment number of events
          nb_events++;
        }
        // Increment time block index
        time_block_index++;
      }
    }

    // Verbose
    if (verbose>1) cout << "--> Number of events: " << nb_events << endl;

    // Close castor datafile
    if (flag_write)
    {
      fclose(fcastor);
      if (nb_data_written != nb_events * 3)
      {
          cerr << "***** Failed to write all data in castor data file !" << endl;
          return 1;
      }
      if (verbose>1) cout << "--> CASToR datafile '" << castor_data_file << "' written" << endl;
    }
    // Write castor header file
    if (flag_write)
    {
      string castor_header_file = output_base_name + ".cdh";
      ofstream cdh(castor_header_file.c_str());
      if (!cdh)
      {
        cerr << "***** Failed to create castor header file '" << castor_header_file << "' for writing !" << endl;
        exit(1);
      }
      cdh << "Scanner name: " << scanner_name << endl;
      cdh << "Data filename: " << castor_data_file << endl;
      cdh << "Number of events: " << nb_events << endl;
      cdh << "Data mode: listmode" << endl;
      cdh << "Data type: PET" << endl;
      cdh << "Start time (s): 0" << endl;
      cdh << "Duration (s): 1" << endl;
      // TODO: time management
      cdh.close();
      // Verbose
      if (verbose>1) cout << "--> CASToR header file '" << castor_header_file << "' written" << endl;
    }

  // End if going through listmode data
  }

  // End
  return 0;
}
