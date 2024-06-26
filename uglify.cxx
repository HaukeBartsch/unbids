#include <cstddef>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataObject.h"
#include "itkSpatialOrientation.h"

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

#include "itkOrientImageFilter.h"

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
#include "itkExtractImageFilter.h"

#include "itkGDCMImageIO.h"

#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <codecvt>
#include <locale> // wstring_convert
#include <map>
#include "dcmtk/dcmdata/dcuid.h"
#include <filesystem>
namespace fs = std::filesystem;

#include <itkNiftiImageIO.h>

// this was useful in creating this program
// https://stackoverflow.com/questions/39150074/generating-a-simple-ct-dicom-image-using-gdcm
 

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

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t") {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange); 
}

typedef itk::Image<double, 3>  DWI;
typedef itk::Image<double, 4>  DWI4;

template< class TImageType = DWI >
std::pair< std::string, typename TImageType::Pointer > GetImageOrientation(const typename TImageType::Pointer inputImage, const std::string &desiredOrientation = "RAI") {
  if (TImageType::ImageDimension != 3) {
    std::cerr << "This function is only defined for 3D images.\n";
    exit(EXIT_FAILURE);
  }
  auto orientFilter = itk::OrientImageFilter< TImageType, TImageType >::New();
  orientFilter->SetInput(inputImage);
  orientFilter->UseImageDirectionOn();
  orientFilter->SetDirectionTolerance(0);
  orientFilter->SetCoordinateTolerance(0);

  auto desiredOrientation_wrap = desiredOrientation;
  std::transform(desiredOrientation_wrap.begin(), desiredOrientation_wrap.end(), desiredOrientation_wrap.begin(), ::toupper);
  
  std::map< std::string, itk::SpatialOrientationEnums::ValidCoordinateOrientations > orientationMap;
  orientationMap["Axial"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RAI;
  orientationMap["Coronal"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RSA;
  orientationMap["Sagittal"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ASL;
  orientationMap["RIP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RIP;
  orientationMap["LIP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LIP;
  orientationMap["RSP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RSP;
  orientationMap["LSP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LSP;
  orientationMap["RIA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RIA;
  orientationMap["LIA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LIA;
  orientationMap["RSA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RSA;
  orientationMap["LSA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LSA;
  orientationMap["IRP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IRP;
  orientationMap["ILP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ILP;
  orientationMap["SRP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SRP;
  orientationMap["SLP"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SLP;
  orientationMap["IRA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IRA;
  orientationMap["ILA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ILA;
  orientationMap["SRA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SRA;
  orientationMap["SLA"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SLA;
  orientationMap["RPI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RPI;
  orientationMap["LPI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LPI;
  orientationMap["RAI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RAI;
  orientationMap["LAI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LAI;
  orientationMap["RPS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RPS;
  orientationMap["LPS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LPS;
  orientationMap["RAS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_RAS;
  orientationMap["LAS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_LAS;
  orientationMap["PRI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PRI;
  orientationMap["PLI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PLI;
  orientationMap["ARI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ARI;
  orientationMap["ALI"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ALI;
  orientationMap["PRS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PRS;
  orientationMap["PLS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PLS;
  orientationMap["ARS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ARS;
  orientationMap["ALS"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ALS;
  orientationMap["IPR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IPR;
  orientationMap["SPR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SPR;
  orientationMap["IAR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IAR;
  orientationMap["SAR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SAR;
  orientationMap["IPL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IPL;
  orientationMap["SPL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SPL;
  orientationMap["IAL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_IAL;
  orientationMap["SAL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_SAL;
  orientationMap["PIR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PIR;
  orientationMap["PSR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PSR;
  orientationMap["AIR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_AIR;
  orientationMap["ASR"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ASR;
  orientationMap["PIL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PIL;
  orientationMap["PSL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_PSL;
  orientationMap["AIL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_AIL;
  orientationMap["ASL"] = itk::SpatialOrientationEnums::ValidCoordinateOrientations::ITK_COORDINATE_ORIENTATION_ASL;

  // set the desired orientation and update
  orientFilter->SetDesiredCoordinateOrientation(orientationMap[desiredOrientation_wrap]);
  orientFilter->Update();
  auto outputImage = orientFilter->GetOutput();

  std::string originalOrientation;

  for (auto it = orientationMap.begin(); it != orientationMap.end(); ++it) {
    if (it->second == orientFilter->GetGivenCoordinateOrientation()) {
      originalOrientation = it->first;
    }
  }
  if (originalOrientation.empty()) {
    originalOrientation = "Unknown";
  }

  return std::make_pair(originalOrientation, outputImage);
}

void convert(json data, std::string nifti_file, std::string output_folder, std::string identifier, std::string StudyInstanceUID, std::string frameOfReferenceUID, int SeriesNumber, bool isMask, bool verbose) {
  // parse the json structure
  //for (auto& [key, value] : data.items()) {
  //  std::cout << key << " : " << value << "\n";
  //}

  // TODO: add an -u option for static uids based on the input folder (md5 of the nii.gz?)
  gdcm::UIDGenerator fuid;
  fuid.SetRoot("1.3.6.1.4.1.45037.1"); // we add a post-fix to the root to indicate its not created by research PACS
  std::string SeriesInstanceUID = fuid.Generate();
  if (verbose)
    fprintf(stdout, "\t  StudyInstanceUID: %s\n\t  SeriesInstanceUID: %s [%d]\n\t  %s\t", StudyInstanceUID.c_str(), SeriesInstanceUID.c_str(), SeriesNumber, output_folder.c_str());

  // read the image data from the nii.gz or .nii file
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(nifti_file.c_str(),  itk::CommonEnums::IOFileMode::ReadMode);

  imageIO->SetFileName(nifti_file);
  imageIO->ReadImageInformation(); // this would work for non .nii files as well. We should support all filetypes (nrrd)
  //itk::CommonEnums::IOPixel pixel_type = imageIO->GetPixelType();
  int dims = (int)imageIO->GetNumberOfDimensions();

  itk::CommonEnums::IOComponent ii= imageIO->GetComponentType();

  time_t studydatetime = gdcm::System::FileTime( nifti_file.c_str() );
  char date[22];
  gdcm::System::FormatDateTime(date, studydatetime);

  typedef itk::ImageFileReader<DWI> DWIReader;
  std::vector<DWI::Pointer> volumes; // for a single series, assume time
  // We might get 4D data here as well. In that case we need to extract individual frames (x,y,z,channel)
  if (dims == 3 /* && imageIO->GetComponentType() == imageIO->GetComponentTypeFromString("double") */) {
    DWIReader::Pointer dwi_reader = DWIReader::New();
    dwi_reader->SetFileName(imageIO->GetFileName());
    dwi_reader->Update();

    /*DWI::RegionType region;
    DWI::IndexType start;
    DWI::SizeType size;

    for (int i = 0; i < dims; ++i){
      start[i] = 0;
      size[i] = imageIO->GetDimensions(i);
    }

    region.SetSize(size);
    region.SetIndex(start); */
    // do we need that now?
    std::pair< std::string, typename DWI::Pointer > erg = GetImageOrientation(dwi_reader->GetOutput());

    DWI::Pointer dwi = erg.second; // dwi_reader->GetOutput();
    volumes.push_back(dwi);
  } else if (dims == 4) {
    // so its 4D? the brats data is using X, Y, Z, N (with N = 4)
    // Can we assume that its always that dimension? Can we assume 
    // that the smallest N is t?
    // Extract using the ExtractImageFilter
    typedef itk::ImageFileReader<DWI4> DWI4Reader;

    DWI4Reader::Pointer dwi_reader = DWI4Reader::New();
    dwi_reader->SetFileName(imageIO->GetFileName());
    dwi_reader->Update();
    DWI4::RegionType inputRegion = dwi_reader->GetOutput()->GetLargestPossibleRegion();

    using ExtractFilterType = itk::ExtractImageFilter<DWI4, DWI>;
    auto extractFilter = ExtractFilterType::New();

    extractFilter->SetDirectionCollapseToSubmatrix();
    DWI4::SizeType roi_size = inputRegion.GetSize();
    int numTs = roi_size[3]; // assume its the 4th dimension
    for (int t = 0; t < numTs; t++) {
      roi_size[3] = 0; // collapse this dimension
      DWI4::IndexType roi_start = inputRegion.GetIndex();
      roi_start[0] = 0;
      roi_start[1] = 0;
      roi_start[2] = 0;
      roi_start[3] = t;

      DWI4::RegionType desiredRegion;
      desiredRegion.SetSize(roi_size);
      desiredRegion.SetIndex(roi_start);
      extractFilter->SetInput(dwi_reader->GetOutput());
      extractFilter->SetExtractionRegion(desiredRegion);
      extractFilter->Update();
      DWI::Pointer p = extractFilter->GetOutput();
      p->DisconnectPipeline(); // remove from the pipeline to get another volume next iteration
      std::pair< std::string, typename DWI::Pointer > erg = GetImageOrientation(p);
      DWI::Pointer dwi = erg.second; // dwi_reader->GetOutput();
      volumes.push_back( p );
    }
  } else {
    fprintf(stderr, "Error: NIfTI is not 3D but %dD.\n", dims);
  }

  //gdcm::UIDGenerator fuid;
  //fuid.SetRoot("1.3.6.1.4.1.45037");
  //std::string SeriesInstanceUID = fuid.Generate();

  for (int vol = 0; vol < volumes.size(); vol++) {
    DWI::Pointer dwi = volumes[vol];
    DWI::RegionType region = dwi->GetLargestPossibleRegion();
    DWI::SizeType size = region.GetSize();

    // if we are not using a mask we should scale the output to a good range
    float minValue = 0.0;
    float maxValue = 1.0;

    // compute good min/max values based on cummulative histograms? or use the whole range?
    // probably best to use the whole range and scale intercept correctly so no data gets lost
    itk::ImageRegionIterator<DWI> inputIterator3D(dwi, dwi->GetLargestPossibleRegion());
    inputIterator3D.GoToBegin();
    bool startD = true;
    while (!inputIterator3D.IsAtEnd()) {
      float val = inputIterator3D.Get();
      if (startD) {
        minValue = (float)val;
        maxValue = (float)val;
        startD = false;
      } else {
        if (val < minValue)
          minValue = (float)val;
        if (val > maxValue)
          maxValue = (float)val;
      }
      ++inputIterator3D;
    }
    if (verbose)
      fprintf(stdout, "[%lu files, %f..%f]\n", size[2], minValue, maxValue);


    // 
    // create a series of 2D images and save in output
    //
    // example from https://itk.org/Doxygen50/html/WikiExamples_2DICOM_2ResampleDICOM_8cxx-example.html
    //

    using ImageType = itk::Image< unsigned short, 2 >;
    auto InputImagePositionPatient = json::array();
    if (data.contains("ImagePositionPatient")) {
       InputImagePositionPatient.push_back(data["ImagePositionPatient"][0]);
       InputImagePositionPatient.push_back(data["ImagePositionPatient"][1]);
       InputImagePositionPatient.push_back(data["ImagePositionPatient"][2]);
    } else {
       InputImagePositionPatient.push_back(0);
       InputImagePositionPatient.push_back(0);
       InputImagePositionPatient.push_back(0);
    }
    for (int f = 0; f < size[2]; f++) { // for each slice f, for each volume vol
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
      if (!data.contains("ImageOrientationPatient") || data["ImageOrientationPatient"].size() != 6) {
        data["ImageOrientationPatient"] = json::array();
        data["ImageOrientationPatient"].push_back(1);
        data["ImageOrientationPatient"].push_back(0);
        data["ImageOrientationPatient"].push_back(0);
        data["ImageOrientationPatient"].push_back(0);
        data["ImageOrientationPatient"].push_back(1);
        data["ImageOrientationPatient"].push_back(0);
      }

      double dirCos[6]={data["ImageOrientationPatient"][0],data["ImageOrientationPatient"][1],data["ImageOrientationPatient"][2],data["ImageOrientationPatient"][3],data["ImageOrientationPatient"][4],data["ImageOrientationPatient"][5]};
      im->SetDirectionCosines( dirCos );
      im->SetSpacing(0, dwi->GetSpacing()[0]);
      im->SetSpacing(1, dwi->GetSpacing()[1]);
      //double pos[3] = {1, 2, 3};
      //im->SetOrigin(0,pos);

      //im->SetTransferSyntax(gdcm::TransferSyntax::ExplicitVRLittleEndian);
      unsigned short *buffer = new unsigned short[size2d[0] * size2d[1] * 2];
      int ii = 0;
      float tmp = (-(minValue*4096.0)/maxValue) / (1.0 - (minValue/maxValue));
      float slope = 1.0/((4096.0 - tmp)/(maxValue));
      float intercept = minValue;
      if (isMask) {
        intercept = 0;
        slope = 1.0;
      }
      // what to do if the slope is too small?
      if (slope < 0.001) {
        slope *= 1000;
      }
      if (!data.contains("RescaleIntercept")) {
        // Intercept and Slope do not work if modality is OT, switch to a better modality like MR
        if (!data.contains("Modality")) {
          if (verbose)
            fprintf(stdout, "Warning: no modality defined. Set Modality=MR to be able to store slope/intercept.\n");
          data["Modality"] = "MR";
        }
        data["RescaleIntercept"] = intercept;
      }
      if (!data.contains("RescaleSlope"))
        data["RescaleSlope"] = slope;

      if (!data.contains("SmallestImagePixelValue")) {
        data["SmallestImagePixelValue"] = minValue;
      }
      if (!data.contains("LargestImagePixelValue")) {
        data["LargestImagePixelValue"] = maxValue;
      }
      if (!data.contains("WindowWidth")) {
        // if we use a larger slope our window values are not ok
        float a = 0 * slope + intercept;
        float b = 4096.0 * slope + intercept;
        if (isMask)
          b = 1.0 * slope + intercept;
        data["WindowWidth"] = (b-a);
      }
      if (!data.contains("WindowCenter")) {
        float a = 0 * slope + intercept;
        float b = 4096.0 * slope + intercept;
        if (isMask)
          b = 1.0 * slope + intercept;
        data["WindowCenter"] = a + (b-a)/2.0;
      }

      inputIterator.GoToBegin();
      outputIterator.GoToBegin();
      while (!inputIterator.IsAtEnd()) {
        float val = inputIterator.Get();
        // here we should scale the data to min/max of unsigned short and set the intercept and slope values y = mx + b
        if (!isMask)
          buffer[ii++] = (unsigned short)( (val-minValue)/(maxValue-minValue) * 4096.0 );
        else
          buffer[ii++] = (unsigned short)( val ); // 0 or 1

        ++inputIterator;
      }
      gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
      pixeldata.SetByteValue((char *)buffer, size2d[0] * size2d[1] * 2);
      im->SetDataElement(pixeldata);

      // map 0 4096 to minValue/maxValue

      gdcm::DataElement de3;
      std::ostringstream value;
      std::string val("");

      // lets switch the SOP Class UID based on the modality
      // might not be needed because we are setting modality
/*      if (data.contains("Modality") && data["Modality"] == "CT") {
        const std::string SOP_CLASS_UID = "0008|0016";
        const std::string C_UID = "1.2.840.10008.5.1.4.1.1.2";
        itk::EncapsulateMetaData<std::string>(dict, SOP_CLASS_UID, C_UID);
      } else if (data.contains("Modality") && data["Modality"] == "US") {
        const std::string SOP_CLASS_UID = "0008|0016";
        const std::string C_UID = "1.2.840.10008.5.1.4.1.1.3.1";
        itk::EncapsulateMetaData<std::string>(dict, SOP_CLASS_UID, C_UID);
      } else if (data.contains("Modality") && data["Modality"] == "MR") {
        const std::string SOP_CLASS_UID = "0008|0016";
        const std::string C_UID = "1.2.840.10008.5.1.4.1.1.4";
        itk::EncapsulateMetaData<std::string>(dict, SOP_CLASS_UID, C_UID);
      } else {
        fprintf(stderr, "unknown modality found, keep default SOP Class UID.\n");
      } */

      // direction for slices based on ImageOrientationPatient
      if (data.contains("ImageOrientationPatient")) {
        float offset_dir[3];
        float a1 = (float)data["ImageOrientationPatient"][0];
        float a2 = (float)data["ImageOrientationPatient"][1];
        float a3 = (float)data["ImageOrientationPatient"][2];
        float b1 = (float)data["ImageOrientationPatient"][3];
        float b2 = (float)data["ImageOrientationPatient"][4];
        float b3 = (float)data["ImageOrientationPatient"][5];

        // cross-product = normal vector
        offset_dir[0] = (a2 * b3) - (a3 * b2);
        offset_dir[1] = (a3 * b1) - (a1 * b3);
        offset_dir[2] = (a1 * b2) - (a2 * b1);
        float len = sqrt((offset_dir[0]*offset_dir[0]) + (offset_dir[1]*offset_dir[1]) + (offset_dir[2]*offset_dir[2]));
        offset_dir[0] /= len;
        offset_dir[1] /= len;
        offset_dir[2] /= len;

        // set ImagePositionPatient
        if (!data.contains("SpacingBetweenSlices"))
          data["SpacingBetweenSlices"] = 1.0;
        float v1 = ((float)(InputImagePositionPatient[0]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[0]));
        float v2 = ((float)(InputImagePositionPatient[1]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[1]));
        float v3 = ((float)(InputImagePositionPatient[2]) + (f*(float)(data["SpacingBetweenSlices"]) * offset_dir[2]));
        //data["ImagePositionPatient"] = json::array();
        //data["ImagePositionPatient"].push_back(v1);
        //data["ImagePositionPatient"].push_back(v2);
        //data["ImagePositionPatient"].push_back(v3);
        // ATTENTION: We set the ImagePositionPatient on gdcm::Image
        // ATTENTION: We set the ImageOrientationPatient on gdcm::Image (further up)
        im->SetOrigin(0, v1);
        im->SetOrigin(1, v2);
        im->SetOrigin(2, v3);
        

        //value.str("");
        //value << v1 << "\\" << v2 << "\\" << v3;
        //itk::EncapsulateMetaData<std::string>(dict,"0020|0032", value.str(), );
        //de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x0032));
        //val = zero_pad(value.str());
        //de3.SetByteValue(val.c_str(), val.size());
        //ds.Insert(de3);
      }
      if (!data.contains("SpacingBetweenSlices")) {
        data["SpacingBetweenSlices"] = 1.0;
      }
      data["SliceLocation"] = (float)f * (float)data["SpacingBetweenSlices"];

      //if (data.contains("SpacingBetweenSlices")) {
      //  value.str("");
      //  value << (f * (float)data["SpacingBetweenSlices"]);
      //  // itk::EncapsulateMetaData<std::string>(*dict,"0020|1041", value.str());
      //  de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x1041));
      //  val = zero_pad(value.str());
      //  de3.SetByteValue(val.c_str(), val.size());
      //  ds.Insert(de3);
      //}

      // PatientName
      //de3 = gdcm::DataElement(gdcm::Tag(0x0010,0x0010));
      //val = zero_pad(PatientID);
      //de3.SetByteValue(val.c_str(), val.size());
      //ds.Insert(de3);

      // PatientID
      //de3 = gdcm::DataElement(gdcm::Tag(0x0010,0x0020));
      //val = zero_pad(PatientID);
      //de3.SetByteValue(val.c_str(), val.size());
      //ds.Insert(de3);

      // ReferringPhysician
      //de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0090));
      //val = zero_pad(std::string("EventName:") + EventName);
      //de3.SetByteValue(val.c_str(), val.size());
      //ds.Insert(de3);

      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0023)); // ContentDate
      de3.SetByteValue(date, 8);
      ds.Insert(de3);

      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0033)); // ContentTime
      de3.SetByteValue(date+8, 6);
      ds.Insert(de3);

      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0020)); // StudyDate
      de3.SetByteValue(date, 8);
      ds.Insert(de3);

      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0030)); // StudyTime
      de3.SetByteValue(date+8, 6);
      ds.Insert(de3);

      //itk::EncapsulateMetaData<std::string>(dict,"0020|000d", StudyInstanceUID);
      de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x000d));
      val = StudyInstanceUID;
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      //itk::EncapsulateMetaData<std::string>(dict,"0020|000e", SeriesInstanceUID);
      de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x000e));
      val = SeriesInstanceUID;
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      //itk::EncapsulateMetaData<std::string>(dict,"0020|0052", frameOfReferenceUID);
      de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x0052));
      val = frameOfReferenceUID;
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      std::string sopInstanceUID = fuid.Generate();
      //itk::EncapsulateMetaData<std::string>(dict,"0008|0018", sopInstanceUID);
      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x0018));
      val = sopInstanceUID;
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      //itk::EncapsulateMetaData<std::string>(dict,"0002|0003", sopInstanceUID);
      /*de3 = gdcm::DataElement(gdcm::Tag(0x0002,0x0003));
      val = zero_pad(sopInstanceUID);
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3); */

      value.str("");
      value << (size[2] * vol) + f + 1; // InstanceNumber
      //itk::EncapsulateMetaData<std::string>(dict,"0020|0013", value.str());
      de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x0013)); // InstanceNumber
      val = zero_pad(value.str());
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

/*
      de3 = gdcm::DataElement(gdcm::Tag(0x0028,0x1052)); // RescaleIntercept
      val = zero_pad(std::to_string(intercept));
      de3.SetByteValue(val.c_str(), val.size());
      de3.SetVR(gdcm::VR::DS);
      ds.Insert(de3);


      //value.str("");
      //value << std::to_string(intercept);
      //itk::EncapsulateMetaData<std::string>(dict,"0028|1052", value.str());

      value.str("");
      value << std::to_string(slope);
      //itk::EncapsulateMetaData<std::string>(dict,"0028|1053", value.str());
      de3 = gdcm::DataElement(gdcm::Tag(0x0028,0x1053)); // RescaleSlope
      val = zero_pad(value.str());
      de3.SetByteValue(val.c_str(), val.size());
      de3.SetVR(gdcm::VR::DS);
      ds.Insert(de3);
*/
      // Series Description - Append new description to current series
      // description
      std::string oldSeriesDesc = "";
      value.str("");
      value << oldSeriesDesc << std::string("uglified ") << identifier;
      // This is an long string and there is a 64 character limit in the
      // standard
      unsigned lengthDesc = value.str().length();
      std::string seriesDesc( value.str(), 0,
                              lengthDesc > 64 ? 64
                              : lengthDesc);
      //itk::EncapsulateMetaData<std::string>(dict,"0008|103e", seriesDesc);
      de3 = gdcm::DataElement(gdcm::Tag(0x0008,0x103e));
      val = zero_pad(seriesDesc);
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      // Series Number
      value.str("");
      value << SeriesNumber;
      //itk::EncapsulateMetaData<std::string>(dict,"0020|0011", value.str());
      de3 = gdcm::DataElement(gdcm::Tag(0x0020,0x0011));
      val = zero_pad(value.str());
      de3.SetByteValue(val.c_str(), val.size());
      ds.Insert(de3);

      const gdcm::Global& g = gdcm::Global::GetInstance();
      const gdcm::Dicts &dicts = g.GetDicts();
      const gdcm::Dict &pubdict = dicts.GetPublicDict();

      // parse the json data and add as new tags before writing
      for (auto& [key, value] : data.items()) {

        if (key == "ImagePositionPatient" || key == "ImageOrientationPatient") // ignore these as they are set on gdcm::Image (im->SetDirectionCosines, im->SetOrigin)
          continue; // ignore, its already set above

        gdcm::Tag t;
        gdcm::DictEntry ent = pubdict.GetDictEntryByKeyword(key.c_str(), t); // slow
        if (ent.GetName() == std::string(""))
          continue;
        gdcm::DataElement de = gdcm::DataElement(gdcm::Tag(t.GetGroup(),t.GetElement()));

        // test if value is a float or string
        if (value.type() == nlohmann::detail::value_t::string) {
          std::string val = zero_pad(std::string(value));
          // de.SetVR( gdcm::Tag(t.GetGroup(),t.GetElement()).GetVR() );
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          de.SetVR( entry.GetVR() );
          //if (key == "PatientPosition") { // this is CS VR, does not work if we don't set the right VR on the de
          //  de.SetVR( gdcm::Attribute<0x0018,0x5100>::GetVR() );
          //}
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::number_float) {
          // this is likely a VR::FD so we need to create such a field
          // see https://github.com/malaterre/GDCM/blob/master/Source/MediaStorageAndFileFormat/gdcmJSON.cxx
          //std::string val = std::to_string((float)value);
          //val = zero_pad(val);
          //de.SetByteValue(val.c_str(), val.size());
          gdcm::DataElement locde;
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          if (entry.GetVR() == gdcm::VR::FD) { // doubles
            gdcm::Element<gdcm::VR::FD,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( 1 * vrsizeof );
            const double v = (double)value;
            el.SetValue(v, 0);
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else if (entry.GetVR() == gdcm::VR::DS && entry.GetVM() == gdcm::VM::VM1_n) {
            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( 1 * vrsizeof );
            const double v = (double)value;
            el.SetValue(v);
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else if (entry.GetVR() == gdcm::VR::DS && entry.GetVM() == gdcm::VM::VM1) {
            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM1> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            //el.SetLength( 1 * vrsizeof );
            const double v = (double)value;
            el.SetValue(v);
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          }

        } else if (value.type() == nlohmann:: detail::value_t::number_integer) {
          //std::string val = std::to_string((int)value);
          //val = zero_pad(val);
          //de.SetByteValue(val.c_str(), val.size());
          gdcm::DataElement locde;
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          if (entry.GetVR() == gdcm::VR::US_SS) { // unsigned short
            gdcm::Element<gdcm::VR::US,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( 1 * vrsizeof );
            const double v = (double)value;
            el.SetValue(v, 0);
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( gdcm::VR::US );
            }
          }
        } else if (value.type() == nlohmann:: detail::value_t::number_unsigned) {
          std::string val = std::to_string((unsigned int)value);
          val = zero_pad(val);
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          de.SetVR( entry.GetVR() );
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann:: detail::value_t::boolean) {
          std::string val = std::to_string((bool)value);
          val = zero_pad(val);
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          de.SetVR( entry.GetVR() );
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_number_float(); })) {
          std::vector<float> ar = value;
          gdcm::DataElement locde;
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          if (entry.GetVR() == gdcm::VR::FD && entry.GetVM() == gdcm::VM::VM1_n) { // doubles
            //gdcm::Element el = entry.GetElement();
            gdcm::Element<gdcm::VR::FD,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const double v = (double)ar[i];
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else if (entry.GetVR() == gdcm::VR::DS && entry.GetVM() == gdcm::VM::VM6) { // doubles
            //gdcm::Element el = entry.GetElement();
            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM6> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            //el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const double v = (double)ar[i];
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
              //de.SetVM( entry.GetVM() );
            }
          } else if (entry.GetVR() == gdcm::VR::FD && entry.GetVM() == gdcm::VM::VM3) { // doubles
            //gdcm::Element el = entry.GetElement();
            gdcm::Element<gdcm::VR::FD,gdcm::VM::VM3> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            //el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const double v = (double)ar[i];
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else if (entry.GetVR() == gdcm::VR::DS && entry.GetVM() == gdcm::VM::VM3) { // doubles
            //gdcm::Element el = entry.GetElement();
            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM3> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            //el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const double v = (double)ar[i];
              el.SetValue((float)ar[i], i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetVLToUndefined();
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else if (entry.GetVR() == gdcm::VR::DS) {
            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const double v = (double)ar[i];
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          } else {
            fprintf(stdout, "Error: cannto write this type.\n");
          }

          //std::string val = "";
          //for (int i = 0; i < ar.size(); i++) {
          //  val += std::to_string(ar[i]);
          //  if (i < ar.size()-1)
          //    val += "\\\\";
          //}
          //val = zero_pad(val);
          //de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_number_integer(); })) {
          std::vector<int> ar = value;
          gdcm::DataElement locde;
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          if (entry.GetVR() == gdcm::VR::US) { // unsigned int as in AcquisitionMatrix
            gdcm::Element<gdcm::VR::US,gdcm::VM::VM1_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const unsigned int v = (unsigned int)ar[i];
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          }
          //std::string val = "";
          //for (int i = 0; i < ar.size(); i++) {
          //  val += std::to_string(ar[i]);
          //  if (i < ar.size()-1)
          //    val += "\\\\";
          //}
          //val = zero_pad(val);
          //de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_number_unsigned(); })) {
          std::vector<unsigned int> ar = value;
          std::string val = "";
          for (int i = 0; i < ar.size(); i++) {
            val += std::to_string(ar[i]);
            if (i < ar.size()-1)
              val += "\\";
          }
          val = zero_pad(val);
          de.SetByteValue(val.c_str(), val.size());
        } else if (value.type() == nlohmann::detail::value_t::array && 
                 std::all_of(value.begin(), value.end(), [](const json& el){ return el.is_string(); })) {
          std::vector<std::string> ar = value;

          gdcm::DataElement locde;
          const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
          if (entry.GetVR() == gdcm::VR::CS) { // unsigned int as in AcquisitionMatrix
            gdcm::Element<gdcm::VR::CS,gdcm::VM::VM2_n> el;
            const int vrsizeof = (entry.GetVR() == gdcm::VR::INVALID ? 0 : entry.GetVR().GetSizeof());
            el.SetLength( ar.size() * vrsizeof );
            for (int i = 0; i < ar.size(); i++) {
              const std::string v = trim(ar[i]);
              //v = std::regex_replace(v, std::regex("^ +| +$|( ) +"), "$1");
              el.SetValue(v, i);
            }
            locde = el.GetAsDataElement();
            if (!locde.IsEmpty()) {
              de.SetValue( locde.GetValue() );
              de.SetVR( entry.GetVR() );
            }
          }

          //std::string val = "";
          //for (int i = 0; i < ar.size(); i++) {
          //  val += std::string(ar[i]);
          //  if (i < ar.size()-1)
          //    val += "\\";
          //}
          //val = zero_pad(val);
          //de.SetByteValue(val.c_str(), val.size());
        } else { // assume a string
          fprintf(stderr, "Error: unsupported type encountered, ignore key: %s\n", key.c_str());
        }
        ds.Insert(de);
      }

      // now save this image slice to a file in the output directory
      //using WriterType = itk::ImageFileWriter<ImageType>;
      std::string ts("");
      if (volumes.size() > 1) {
        ts = "_" + leading_zeros(std::to_string(vol),2);
      }
      std::string output_fname = output_folder + std::string("/image_") + identifier + ts + std::string("_") + leading_zeros(std::to_string(f),4) + std::string(".dcm");

      gdcm::ImageWriter w;
      w.SetImage(*im);
      w.SetFile(*filePtr);
      w.SetFileName(output_fname.c_str());
      if (!w.Write()) {
        fprintf(stderr, "ERROR: writing file\n");
      }
    }
  } // for volumes

}

std::string randomHexInt(int N) { 
    srand(time(0)); 

    int maxSize = 2; 
    // Stores all the possible characters 
    // in the Hexadecimal notation 
    char hexChar[] 
        = { '0', '1', '2', '3', '4', '5', 
            '6', '7', '8', '9', 'A', 'B', 
            'C', 'D', 'E', 'F' }; 
    std::ostringstream value;
    value.str("");

    // Loop to print N integers 
    for (int i = 0; i < N; i++) { 
  
        // Randomly select length of the 
        // int in the range [1, maxSize] 
        int len = rand() % maxSize + 1; 
        // first character should not be 0

        // Print len characters 
        for (int j = 0; j < len; j++) { 
  
            // Print a randomly selected 
            // character 
            value << hexChar[rand() % 16]; 
        }  
    } 
    return value.str();
} 

// We would like to call this with a single case. The raw data should come from -i, apply to all files in that folder.
// The mask data should come from another folder (derivatives).
// Usage: 
//      -i sub-strokecase0001 -m sub-strokecase0001
// Each of the -i nii.gz files will have a corresponding .json. The -m nii.gz file shall borrow its json
// from a matching volume of -i. 

// TODO: add a -m option to modifiy DICOM tags. This does only make sense if 
// we can provide it several times (GetValueAsLIST).

namespace po = boost::program_options;



int main(int argc, char *argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");

  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();
  resultJSON["run_date_time"] = to_simple_string(timeLocal);

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  bool verbose = false;
  std::string output("");
  json overwrites; // command line options to overwrite entries derived from side-loading json

  po::options_description desc("UGLIFY: Convert a volume file to a DICOM folder.\nExample:\n  uglify --verbose -o /tmp/b4 --modify PatientID=MEME data/brats/BRATS_280.nii.gz\nAllowed options");
  desc.add_options()
      ("help,h", "Print this help.")
      ("version,V", "Print the version number.")
      ("verbose,v", po::bool_switch(&verbose), "Print more verbose output during processing.")
      ("raw-data,i", po::value< std::vector<std::string> >(), "Folder with nii.gz and .json files, or a single .nii/.nii.gz file.")
      ("mask-data,m", po::value< std::vector<std::string> >(), "Folder with nii.gz and .json files, or a single .nii/.nii.gz file.")
      ("outdir,o", po::value<std::string>(&output), "Output directory for DICOM files")
      ("modify,e", po::value< std::vector< std::string > >(), "Modify individual DICOM tags (\"PatientID=MEME\"). This option can be used more than once.")
  ;
  // allow positional arguments to map to rawdata
  po::positional_options_description p;
  p.add("raw-data", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);
  } catch(std::exception& e) {
        std::cout << e.what() << "\n";
        return 1;
  }

  if (vm.count("help") || argc <= 1) {
    std::cout << desc << std::endl;
    return 0;
  }

  std::string versionString = std::string("0.0.1.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  if (vm.count("version")) {
    fprintf(stdout, "version: %s\n", versionString.c_str());
    return 0;
  }

  if (!vm.count("outdir")) {
    fprintf(stderr, "Error: we require an output directory -o.\n");
    return -1;
  }

  std::vector< std::string > image_files;
  std::vector<std::string> image_folders;
  if (vm.count("raw-data")) {
    auto rd = vm["raw-data"].as< std::vector< std::string > >();
    for (int i = 0; i < rd.size(); i++) {
      std::string tmp = rd[i];
      if (itksys::SystemTools::FileIsDirectory(tmp.c_str())) {
        image_folders.push_back(tmp);
      } else if (std::filesystem::is_regular_file(tmp.c_str())) {
        image_files.push_back(tmp);
      } else {
        fprintf(stderr, "Error: raw data argument is not a valid directory or file [%s]\n", tmp.c_str());
        exit(-1);
      }
    }
  }

  std::vector< std::string > mask_files;
  std::vector<std::string> mask_folders;
  if (vm.count("mask-data")) {
    auto rd = vm["mask-data"].as< std::vector< std::string > >();
    for (int i = 0; i < rd.size(); i++) {
      std::string tmp = rd[i];
      if (itksys::SystemTools::FileIsDirectory(tmp.c_str())) {
        mask_folders.push_back(tmp);
      } else if (std::filesystem::is_regular_file(tmp.c_str())) {
        mask_files.push_back(tmp);
      } else {
        fprintf(stderr, "Error: raw data argument is not a valid directory or file [%s]\n", tmp.c_str());
        exit(-1);
      }
    }
  }

  if (vm.count("modify")) {
    auto md = vm["modify"].as< std::vector< std::string > >();
    for (int i = 0; i < md.size(); i++) {
      //fprintf(stdout, "Found a modify option: %s\n", md[i].c_str());
      // we should add these to the data json
      std::vector<std::string> pieces;
      tokenize(md[i], '=', pieces);
      if (pieces.size() == 2) {
        overwrites[pieces[0]] = pieces[1];
      } else {
        fprintf(stderr, "Error: could not parse modify option \"%s\" as \"key=value\".", md[i].c_str());
        exit(-1);
      }
    }
  }

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }

  gdcm::UIDGenerator fuid;
  fuid.SetRoot("1.3.6.1.4.1.45037.1"); // add .1 to our root to mark files as non-rPACS
  std::string StudyInstanceUID = fuid.Generate();
  std::string frameOfReferenceUID = fuid.Generate();
  int SeriesCounter = 1;
  std::string json_dummy_file = "";
  // create an AccessionNumber and a StudyID (16 characters each)
  // they are shared for all files in this study
  std::string AccessionNumber = randomHexInt(8);
  std::string StudyID = randomHexInt(8);
  std::string PatientID = randomHexInt(8);
  std::string StudyDate = "";
  std::string StudyTime = "";
  std::string EventName("unknown");

  // we should start by parsing the image_folders for nii.gz files (and the corresponding .json)
  for (int i = 0; i < image_folders.size(); i++) {
    std::string path = image_folders[i];
    for (const auto & entry: fs::recursive_directory_iterator(path)) {
      std::string fn = entry.path().string();
      if (std::filesystem::is_regular_file(fn)) { // ignore directories
        // ignore all files that are not .nii or .nii.gz
        if (!ends_with(fn, ".nii.gz") && !ends_with(fn, ".nii"))
          continue; // ignore
        image_files.push_back(fn);
      }
    }
  }

  json json_data; // do we have to reset this every time? (no, we are reading from disk)
  for (int i = 0; i < image_files.size(); i++) {
    json_data = json::object();

    std::string fn = image_files[i];
    // if we have a .nii.gz or .nii we can start
    std::string identifier("unknown");
    std::string json_file = fn;
    fs::directory_entry entry(fn);
    std::string extension = "";
    if (ends_with(fn, ".nii.gz")) {
      std::string f_only = entry.path().filename().string();
      std::vector<std::string> pieces;
      tokenize(f_only.substr(0, f_only.size() - std::string(".nii.gz").size()), '_', pieces);
      if (pieces.size() > 2) {
        identifier = pieces[2]; // something like "adc"
      }
      if (pieces.size() > 0) {
        // remove the sub- component if it exists
        if (pieces[0].substr(0,4) == std::string("sub-")) {
          PatientID = pieces[0].substr(4, pieces[0].size());
        }
      }
      if (pieces.size() > 1) {
        // remove the sub- component if it exists
        if (pieces[1].substr(0,4) == std::string("ses-")) {
          EventName = pieces[1].substr(4, pieces[1].size());
        }
      }
      json_file = json_file.substr(0, json_file.size() - std::string(".nii.gz").size()) + ".json";
      extension = ".gz";
    } else {
      std::string f_only = entry.path().filename().string();
      std::vector<std::string> pieces;
      tokenize(f_only.substr(0, f_only.size() - std::string(".nii").size()), '_', pieces);
      if (pieces.size() > 2) {
        identifier = pieces[2]; // something like "adc"
      }
      if (pieces.size() > 0) {
        // remove the sub- component if it exists
        if (pieces[0].substr(0,4) == std::string("sub-")) {
          PatientID = pieces[0].substr(4, pieces[0].size());
        }
      }
      if (pieces.size() > 1) {
        // remove the sub- component if it exists
        if (pieces[1].substr(0,4) == std::string("ses-")) {
          EventName = pieces[1].substr(4, pieces[1].size());
        }
      }
      json_file = json_file.substr(0, json_file.size() - std::string(".nii").size()) + ".json";
    }
    // check if that the json file exists, if not use a dummy json instead
    if (std::filesystem::is_regular_file(json_file)) {
      json_dummy_file = json_file; // keep a record of this for the masks
      // we found an nii and a corresponding json file
      if (verbose)
        fprintf(stdout, "%s\n%s\n", fn.c_str(), json_file.c_str() );
        // read the json and start processing
      std::ifstream f(json_file);
      json_data = json::parse(f);
      // add more tags to the json_data before writing
      // TODO: check if they exist and only add if not
      json_data["PatientID"] = PatientID;
      json_data["PatientName"] = PatientID;
      json_data["ReferringPhysicianName"] = EventName;
      json_data["AccessionNumber"] = AccessionNumber;
      json_data["StudyID"] = StudyID;
    } else {          
      if (verbose)
        fprintf(stdout, "%s\n", fn.c_str());
      json_data["PatientID"] = PatientID;
      json_data["PatientName"] = PatientID;
      json_data["ReferringPhysicianName"] = EventName;
      json_data["AccessionNumber"] = AccessionNumber;
      json_data["StudyID"] = StudyID;
    }
    // add the identifier as a folder
    std::string foldername = output + "/" + PatientID + "/" + EventName + "/" + std::to_string(SeriesCounter) + "_" + identifier + "/";
    if (!itksys::SystemTools::FileIsDirectory(foldername.c_str())) {
        // create the output directory
        create_directories(foldername);
    }
    // add our overwrites
    if (!overwrites.is_null())
      json_data.update(overwrites);

    convert(json_data, fn, foldername, identifier, StudyInstanceUID, frameOfReferenceUID, SeriesCounter++, false, verbose);
  }

  // repeat the same for any mask
  for (int i = 0; i < mask_folders.size(); i++) {
    std::string path = mask_folders[i];
    for (const auto & entry: fs::recursive_directory_iterator(path)) {
      std::string fn = entry.path().string();
      if (std::filesystem::is_regular_file(fn)) {
        // ignore all files that are not .nii or .nii.gz
        if (!ends_with(fn, ".nii.gz") && !ends_with(fn, ".nii"))
          continue; // ignore
        mask_files.push_back(fn);
      }
    }
  }
  for (int i = 0; i < mask_files.size(); i++) {
    json_data = json::object();
    std::string fn = mask_files[i];
    // if we have a .nii.gz or .nii we can start
    std::string identifier("unknown");
    std::string json_file = fn;
    fs::directory_entry entry(fn);
    std::string extension = "";
    if (ends_with(fn, ".nii.gz")) {
      std::string f_only = entry.path().filename().string();
      std::vector<std::string> pieces;
      tokenize(f_only.substr(0, f_only.size() - std::string(".nii.gz").size()), '_', pieces);
      if (pieces.size() > 2) {
        identifier = pieces[2]; // something like "adc"
      }
      if (pieces.size() > 0) {
        // remove the sub- component if it exists
        if (pieces[0].substr(0,4) == std::string("sub-")) {
          PatientID = pieces[0].substr(4, pieces[0].size());
        }
      }
      if (pieces.size() > 1) {
        // remove the sub- component if it exists
        if (pieces[1].substr(0,4) == std::string("ses-")) {
          EventName = pieces[1].substr(4, pieces[1].size());
        }
      }
      json_file = json_file.substr(0, json_file.size() - std::string(".nii.gz").size()) + ".json";
      extension = ".gz";
    } else {
      std::string f_only = entry.path().filename().string();
      std::vector<std::string> pieces;
      tokenize(f_only.substr(0, f_only.size() - std::string(".nii").size()), '_', pieces);
      if (pieces.size() > 2) {
        identifier = pieces[2]; // something like "adc"
      }
      if (pieces.size() > 0) {
        // remove the sub- component if it exists
        if (pieces[0].substr(0,4) == std::string("sub-")) {
          PatientID = pieces[0].substr(4, pieces[0].size());
        }
      }
      if (pieces.size() > 1) {
        // remove the sub- component if it exists
        if (pieces[1].substr(0,4) == std::string("ses-")) {
          EventName = pieces[1].substr(4, pieces[1].size());
        }
      }
      json_file = json_file.substr(0, json_file.size() - std::string(".nii").size()) + ".json";
    }
    if (!std::filesystem::is_regular_file(json_file) && std::filesystem::is_regular_file(json_dummy_file)) {
      json_file = json_dummy_file;
    }
    //json json_data;
    // check if that the file exists
    if (std::filesystem::is_regular_file(json_file)) {
      // we found an nii and a corresponding json file
      if (verbose)
        fprintf(stdout, "%s\n%s\n", fn.c_str(), json_file.c_str() );
      // read the json and start processing
      std::ifstream f(json_file);
      json_data = json::parse(f);
      json_data["PatientID"] = PatientID;
      json_data["PatientName"] = PatientID;
      json_data["ReferringPhysicianName"] = EventName;
      json_data["AccessionNumber"] = AccessionNumber;
      json_data["StudyID"] = StudyID;
    } else {
      if (verbose)
        fprintf(stdout, "%s\n", fn.c_str() );
      json_data["PatientID"] = PatientID;
      json_data["PatientName"] = PatientID;
      json_data["ReferringPhysicianName"] = EventName;
      json_data["AccessionNumber"] = AccessionNumber;
      json_data["StudyID"] = StudyID;
    }
    // add the identifier as a folder
    std::string foldername = output + "/" + PatientID + "/" + EventName + "/" + std::to_string(SeriesCounter) + "_" + identifier + "/";
    if (!itksys::SystemTools::FileIsDirectory(foldername.c_str())) {
        // create the output directory
        create_directories(foldername);
    }
    // add our overwrites
    if (!overwrites.is_null())
      json_data.update(overwrites);

    // convert as mask
    convert(json_data, fn, foldername, identifier, StudyInstanceUID, frameOfReferenceUID, SeriesCounter++, true, verbose);            
  }

  // check if we have already a tracking file (mapping.json)
  std::ofstream mapFile;
  std::string mapping_file = output + "/" + "mapping.csv";
  if (!itksys::SystemTools::FileExists(mapping_file.c_str(), true)) {
    // create and add a header
    mapFile.open(mapping_file);
    mapFile << "subjectid" << "," << "eventname" << "," << "AccessionNumber" << "," << "StudyID" << "," << "StudyInstanceUID" << std::endl;
    mapFile.close();
  }
  mapFile.open(mapping_file, std::ios_base::app);
  mapFile << json_data["PatientID"] << "," << json_data["ReferringPhysicianName"] << "," << json_data["AccessionNumber"] << "," << json_data["StudyID"] << "," << StudyInstanceUID << std::endl;
  mapFile.close();


  if (verbose) {
    std::string res = resultJSON.dump(4) + "\n";
    fprintf(stdout, "%s", res.c_str());
  }

  return EXIT_SUCCESS;
}
