#include <cstddef>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

#include "itkMinimumMaximumImageCalculator.h"
#include "itkRGBPixel.h"
#include "itkScalarImageToHistogramGenerator.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"

#include "gdcmAnonymizer.h"
#include "gdcmAttribute.h"
#include "gdcmDataSetHelper.h"
#include "gdcmDirectoryHelper.h"
#include "gdcmFileDerivation.h"
#include "gdcmFileExplicitFilter.h"
#include "gdcmGlobal.h"
#include "gdcmImageApplyLookupTable.h"
#include "gdcmImageChangePlanarConfiguration.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmMediaStorage.h"
#include "gdcmReader.h"
#include "gdcmRescaler.h"
#include "gdcmStringFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkConstantPadImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>
#include "dcmtk/dcmdata/dcuid.h"
#include <filesystem>
namespace fs = std::filesystem;

#include <itkNiftiImageIO.h>
 

using json = nlohmann::json;
using namespace boost::filesystem;

json resultJSON;

inline bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline std::string zero_pad(std::string str) {
  if (str.size() % 2 != 0) {
    return str + std::string(" ");
  }
  return str;
}

void tokenize(std::string const &str, const char delim, std::vector<std::string> &out) { 
  size_t start; size_t end = 0;

  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    out.push_back(str.substr(start, end - start));
  }
}

inline std::string leading_zeros(std::string a, int num) {
  return std::string(num - std::min(num, (int)a.length()), '0') + a;
}

void convert(json data, std::string nifti_file, std::string output_folder, std::string StudyInstanceUID, std::string frameOfReferenceUID, int SeriesNumber) {
  // parse the json structure
  //for (auto& [key, value] : data.items()) {
  //  std::cout << key << " : " << value << "\n";
  //}

  // TODO: add an -u option for static uids based on the input folder (md5 of the nii.gz?)
  gdcm::UIDGenerator fuid;
  fuid.SetRoot("1.3.6.1.4.1.45037");
  std::string SeriesInstanceUID = fuid.Generate();

  // read the image data from the nii.gz or .nii file
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(nifti_file.c_str(), itk::ImageIOFactory::ReadMode);

  imageIO->SetFileName(nifti_file);
  imageIO->ReadImageInformation();
  itk::ImageIOBase::IOPixelType pixel_type = imageIO->GetPixelType();
  int dims = (int)imageIO->GetNumberOfDimensions();

  itk::CommonEnums::IOComponent ii= imageIO->GetComponentType();

  // we want to use SpacingBetweenSlices 0.712 with ImageOrientationPatient and ImagePositionPatient
  if (dims == 3 && imageIO->GetComponentType() == imageIO->GetComponentTypeFromString("double")) {
    typedef itk::Image<double, 3>  DWI;
    typedef itk::ImageFileReader<DWI> DWIReader;

    DWIReader::Pointer dwi_reader = DWIReader::New();
    dwi_reader->SetFileName(imageIO->GetFileName());
    dwi_reader->Update();

    DWI::RegionType region;
    DWI::IndexType start;
    DWI::SizeType size;

    for (int i = 0; i < dims; ++i){
      start[i] = 0;
      size[i] = imageIO->GetDimensions(i);
    }

    region.SetSize(size);
    region.SetIndex(start);

    DWI::Pointer dwi = dwi_reader->GetOutput();

    gdcm::UIDGenerator fuid;
    fuid.SetRoot("1.3.6.1.4.1.45037");
    std::string SeriesInstanceUID = fuid.Generate();

    // 
    // create a series of 2D images and save in output
    //
    // example from https://itk.org/Doxygen50/html/WikiExamples_2DICOM_2ResampleDICOM_8cxx-example.html
    //

    using ImageType = itk::Image< unsigned short, 2 >;

    for (int f = 0; f < size[2]; f++) { // for each image
      using ImageIOType = itk::GDCMImageIO;
      auto dicomIO = ImageIOType::New(); // we set dictionary values in this dicomIO and add it in the writer

      itk::MetaDataDictionary &dict = dicomIO->GetMetaDataDictionary();
      ImageType::Pointer image = ImageType::New();

      ImageType::IndexType _start;
      _start.Fill(0);
  
      ImageType::SizeType _size;
      _size[0] = size[0];
      _size[1] = size[1];
  
      ImageType::RegionType _region;
      _region.SetSize(_size);
      _region.SetIndex(_start);
  
      image->SetRegions(_region);
      image->Allocate();

      gdcm::File *filePtr = new gdcm::File;
      gdcm::Anonymizer anon;
      anon.SetFile(*filePtr);
      gdcm::DataSet &ds = filePtr->GetDataSet();
  
      // we can copy over one slice by setting up the region correctly
      DWI::RegionType regionExtract = region;
      DWI::SizeType size2d;
      size2d[0] = _size[0]; size2d[1] = _size[1]; size2d[2] = 1;
      regionExtract.SetSize(size2d);
      DWI::IndexType start2d;
      start2d[0] = 0; start2d[1] = 0; start2d[2] = f;
      regionExtract.SetIndex(start2d);
      // copy image pixel from regionExtract to _region
      itk::ImageRegionIterator<DWI> inputIterator(dwi, regionExtract);
      itk::ImageRegionIterator<ImageType> outputIterator(image, _region);

      gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;
      im->SetNumberOfDimensions(2);
      im->SetDimension(0, size2d[0]);
      im->SetDimension(1, size2d[1]);
      gdcm::PixelFormat pf = gdcm::PixelFormat::UINT16;
      im->SetPixelFormat(pf);
      im->SetPhotometricInterpretation(gdcm::PhotometricInterpretation::MONOCHROME2); // change_image.GetPhotometricInterpretation());
      im->GetPixelFormat().SetSamplesPerPixel(1);
      //im->SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
      unsigned short *buffer = new unsigned short[size2d[0] * size2d[1] * 2];
      int ii = 0;

      inputIterator.GoToBegin();
      outputIterator.GoToBegin();
      while (!inputIterator.IsAtEnd()) {
        float val = inputIterator.Get();
        buffer[ii++] = (unsigned short)val;
        ++inputIterator;
      }
      gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
      pixeldata.SetByteValue((char *)buffer, size2d[0] * size2d[1] * 2);
      im->SetDataElement(pixeldata);

      // direction for slices based on ImageOrientationPatient
      float offset_dir[3];
      float a1 = (float)data["ImageOrientationPatient"][0];
      float a2 = (float)data["ImageOrientationPatient"][1];
      float a3 = (float)data["ImageOrientationPatient"][2];
      float b1 = (float)data["ImageOrientationPatient"][3];
      float b2 = (float)data["ImageOrientationPatient"][4];
      float b3 = (float)data["ImageOrientationPatient"][5];

      offset_dir[0] = (a2 * b3) - (a3 * b2);
      offset_dir[1] = (a3 * b1) - (a1 * b3);
      offset_dir[2] = (a1 * b2) - (a2 * b1);
      float len = sqrt((offset_dir[0]*offset_dir[0]) + (offset_dir[1]*offset_dir[1]) + (offset_dir[2]*offset_dir[2]));
      offset_dir[0] /= len;
      offset_dir[1] /= len;
      offset_dir[2] /= len;

      // set ImagePositionPatient
      std::ostringstream value;
      value.str("");
      value << ((float)(data["ImagePositionPatient"][0]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[0])) << "\\" 
            << ((float)(data["ImagePositionPatient"][1]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[1])) << "\\" 
            << ((float)(data["ImagePositionPatient"][2]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[2]));
      //gdcm::DataElement de2 = gdcm::DataElement(gdcm::Tag(0x0020,0x0032));
      //std::string val = value.str();
      //de2.SetByteValue(val.c_str(), val.size());
      //ds.Insert(de2);

      itk::EncapsulateMetaData<std::string>(dict,"0020|0032", value.str());

      itk::EncapsulateMetaData<std::string>(dict,"0020|000d", StudyInstanceUID);
      itk::EncapsulateMetaData<std::string>(dict,"0020|000e", SeriesInstanceUID);
      itk::EncapsulateMetaData<std::string>(dict,"0020|0052", frameOfReferenceUID);

      std::string sopInstanceUID = fuid.Generate();
      itk::EncapsulateMetaData<std::string>(dict,"0008|0018", sopInstanceUID);
      itk::EncapsulateMetaData<std::string>(dict,"0002|0003", sopInstanceUID);

      value.str("");
      value << f + 1;

      itk::EncapsulateMetaData<std::string>(dict,"0020|0013", value.str());
      // Series Description - Append new description to current series
      // description
      std::string oldSeriesDesc = "";
      value.str("");
      value << oldSeriesDesc << "uglified";
      // This is an long string and there is a 64 character limit in the
      // standard
      unsigned lengthDesc = value.str().length();
      std::string seriesDesc( value.str(), 0,
                              lengthDesc > 64 ? 64
                              : lengthDesc);
      itk::EncapsulateMetaData<std::string>(dict,"0008|103e", seriesDesc);

      // Series Number
      value.str("");
      value << SeriesNumber;
      itk::EncapsulateMetaData<std::string>(dict,"0020|0011", value.str());

      const gdcm::Global& g = gdcm::Global::GetInstance();
      const gdcm::Dicts &dicts = g.GetDicts();
      const gdcm::Dict &pubdict = dicts.GetPublicDict();

      // parse the json data and add as new tags before writing
      for (auto& [key, value] : data.items()) {

        if (key == "ImagePositionPatient")
          continue; // ignore, its already set above

        gdcm::Tag t;
        gdcm::DictEntry ent = pubdict.GetDictEntryByKeyword(key.c_str(), t); // slow
        gdcm::DataElement de = gdcm::DataElement(gdcm::Tag(t.GetGroup(),t.GetElement()));

        // test if value is a float or string
        if (value.type() == nlohmann::detail::value_t::string) {
          std::string val = zero_pad(std::string(value));
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::number_float) {
          std::string val = std::to_string((float)value);
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::number_integer) {
          std::string val = std::to_string((int)value);
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::number_unsigned) {
          std::string val = std::to_string((unsigned int)value);
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::boolean) {
          std::string val = std::to_string((bool)value);
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_number_float(); })) {
          std::vector<float> ar = value;
          std::string val = "";
          for (int i = 0; i < ar.size(); i++) {
            val += std::to_string(ar[i]);
            if (i < ar.size()-1)
              val += "\\\\";
          }
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_number_integer(); })) {
          std::vector<int> ar = value;
          std::string val = "";
          for (int i = 0; i < ar.size(); i++) {
            val += std::to_string(ar[i]);
            if (i < ar.size()-1)
              val += "\\\\";
          }
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_string(); })) {
          std::vector<std::string> ar = value;
          std::string val = "";
          for (int i = 0; i < ar.size(); i++) {
            val += std::string(ar[i]);
            if (i < ar.size()-1)
              val += "\\\\";
          }
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else { // assume a string
          fprintf(stderr, "Error: unsupported type encountered, ignore key: %s\n", key.c_str());
        }
        ds.Insert(de);
      }

      // now save this image slice to a file in the output directory
      //using WriterType = itk::ImageFileWriter<ImageType>;
      std::string output_fname = output_folder + std::string("/image_") + leading_zeros(std::to_string(f),4) + std::string(".dcm");

      gdcm::ImageWriter w;
      w.SetImage(*im);
      w.SetFile(*filePtr);
      w.SetFileName(output_fname.c_str());
      if (!w.Write()) {
        fprintf(stderr, "ERROR: writing file\n");
      }
    }
  }

}

// We would like to call this with a single case. The raw data should come from -i, apply to all files in that folder.
// The mask data should come from another folder (derivatives).
// Usage: 
//      -i sub-strokecase0001 -m sub-strokecase0001
// Each of the -i nii.gz files will have a corresponding .json. The -m nii.gz file shall borrow its json
// from a matching volume of -i. 

int main(int argc, char *argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");

  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  std::string versionString = std::string("0.0.4.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  command.SetVersion(versionString.c_str());
  command.SetDate(to_simple_string(timeLocal).c_str());
  command.SetDescription("UGLIFY: Convert a BIDS like folder to DICOM.");
  command.SetCategory("image conversion");
  command.AddField("outdir", "Directory for output DICOM images.", MetaCommand::STRING, true);

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.SetOptionLongTag("SeriesName", "seriesname");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);

  command.SetOption("RawData", "i", false, "Folder with nii.gz files and .json files.");
  command.SetOptionLongTag("RawData", "raw-data");
  command.AddOptionField("RawData", "value", MetaCommand::STRING, false);

  command.SetOption("MaskData", "m", false, "Folder with nii.gz files representing mask volumes.");
  command.SetOptionLongTag("MaskData", "mask-data");
  command.AddOptionField("MaskData", "value", MetaCommand::STRING, false);


  // allow for interpolation between slices (assumes a single object)
  // instead do this in a separate command (MorphologicalContourInterpolation)
  
  command.SetOption(
      "UIDFixed", "u", false,
      "If enabled identifiers are stable - will not change for a given input. This allows image series to overwrite each other - assuming that the PACS "
      "supports this overwrite mode. By default the SeriesInstanceUID and SOPInstanceUID values are generated again every time the processing is done.");
  command.SetOptionLongTag("UIDFixed", "uid-fixed");

  command.SetOption("Verbose", "v", false, "Print more verbose output");
  command.SetOptionLongTag("Verbose", "verbose");

  command.SetOption("BrightnessContrastLL", "d", false, "Set threshold for brightness / contrast based on cummulative histogram lower limit (percentage dark pixel 0.01).");
  command.SetOptionLongTag("BrightnessContrastLL", "brightness-contrast-ll");
  command.AddOptionField("BrightnessContrastLL", "value", MetaCommand::FLOAT, false);

  command.SetOption("BrightnessContrastUL", "b", false, "Set threshold for brightness / contrast based on cummulative histogram upper limit (percentage bright pixel 0.999).");
  command.SetOptionLongTag("BrightnessContrastUL", "brightness-contrast-ul");
  command.AddOptionField("BrightnessContrastUL", "value", MetaCommand::FLOAT, false);

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  std::vector<std::string> image_folders;
  if (command.GetOptionWasSet("RawData")) {
    image_folders.push_back(command.GetValueAsString("RawData", "value"));
  }

  std::vector<std::string> mask_folders;
  if (command.GetOptionWasSet("MaskData")) {
    mask_folders.push_back(command.GetValueAsString("MaskData", "value"));
  }

  float brightness_contrast_ll = 0.01;
  float brightness_contrast_ul = 0.999;
  float brightnesscontrast_ll = brightness_contrast_ll;
  float brightnesscontrast_ul = brightness_contrast_ul;
  if (command.GetOptionWasSet("BrightnessContrastLL")) {
    brightnesscontrast_ll = command.GetValueAsFloat("BrightnessContrastLL", "value");
    if (brightnesscontrast_ll < 0 || brightnesscontrast_ll > 1.0) {
      fprintf(stdout, "Warning: lower brightness values not between 0 and 1. Adjusted to 0.01.\n");
      brightnesscontrast_ll = 0.01;
    }
  }
  if (command.GetOptionWasSet("BrightnessContrastUL")) {
    brightnesscontrast_ul = command.GetValueAsFloat("BrightnessContrastUL", "value");
    if (brightnesscontrast_ul < 0 || brightnesscontrast_ul > 1.0) {
      fprintf(stdout, "Warning: upper brightness values not between 0 and 1. Adjusted to 0.999.\n");
      brightnesscontrast_ul = 0.999;
    }
  }
  if (brightnesscontrast_ul < brightnesscontrast_ll) {
    float tmp = brightnesscontrast_ll;
    brightnesscontrast_ll = brightnesscontrast_ul;
    brightnesscontrast_ul = tmp;
  }
  brightness_contrast_ll = brightnesscontrast_ll;
  brightness_contrast_ul = brightnesscontrast_ul;
  if (verbose) {
    fprintf(stdout, "create report with brightness/contrast setting %.03f %.03f\n", brightness_contrast_ll, brightness_contrast_ul);
  }

  bool uidFixedFlag = false;
  if (command.GetOptionWasSet("UIDFixed"))
    uidFixedFlag = true;

  bool seriesIdentifierFlag = false;
  std::string output = command.GetValueAsString("outdir");

  if (command.GetOptionWasSet("SeriesName"))
    seriesIdentifierFlag = true;

  std::string seriesName = command.GetValueAsString("SeriesName", "seriesname");

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  gdcm::UIDGenerator fuid;
  fuid.SetRoot("1.3.6.1.4.1.45037");
  std::string StudyInstanceUID = fuid.Generate();
  std::string frameOfReferenceUID = fuid.Generate();
  int SeriesCounter = 1;

  // we should start by parsing the image_folders for nii.gz files (and the corresponding .json)
  for (int i = 0; i < image_folders.size(); i++) {
    std::string path = image_folders[i];
    for (const auto & entry: fs::recursive_directory_iterator(path)) {
      std::string fn = entry.path().string();
      if (std::filesystem::is_regular_file(fn)) {
        // ignore all files that are not .nii or .nii.gz
        if (!ends_with(fn, ".nii.gz") && !ends_with(fn, ".nii"))
          continue; // ignore
        // if we have a .nii.gz or .nii we can start
        std::string identifier("unknown");
        std::string json_file = fn;
        if (ends_with(fn, ".nii.gz")) {
          std::string f_only = entry.path().filename().string();
          std::vector<std::string> pieces;
          tokenize(f_only.substr(0, f_only.size() - std::string(".nii.gz").size()), '_', pieces);
          if (pieces.size() > 2) {
            identifier = pieces[2]; // something like "adc"
          }
          json_file = json_file.substr(0, json_file.size() - std::string(".nii.gz").size()) + ".json";
        } else {
          std::string f_only = entry.path().filename().string();
          std::vector<std::string> pieces;
          tokenize(f_only.substr(0, f_only.size() - std::string(".nii").size()), '_', pieces);
          if (pieces.size() > 2) {
            identifier = pieces[2]; // something like "adc"
          }
          json_file = json_file.substr(0, json_file.size() - std::string(".nii").size()) + ".json";
        }
        // check if that file exists
        if (std::filesystem::is_regular_file(json_file)) {
          // we found an nii and a corresponding json file
          fprintf(stdout, "found an nii file and a matching json:\n\t%s\n\t%s\n", fn.c_str(), json_file.c_str() );

          // read the json and start processing
          std::ifstream f(json_file);
          json json_data = json::parse(f);
          // add the identifier as a folder
          std::string foldername = output + "/" + std::to_string(SeriesCounter) + "_" + identifier + "/";
          if (!itksys::SystemTools::FileIsDirectory(foldername.c_str())) {
              // create the output directory
              create_directories(foldername);
          }

          convert(json_data, fn, foldername, StudyInstanceUID, frameOfReferenceUID, SeriesCounter++);
        }

      }

    }
  }



  if (verbose) {
    std::string res = resultJSON.dump(4) + "\n";
    fprintf(stdout, "%s", res.c_str());
  }

  return EXIT_SUCCESS;
}
