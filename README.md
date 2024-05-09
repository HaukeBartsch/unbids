# un-BIDS

A.I. model training can be done with publicly available data on kaggle etc.. If the data has been converted to nii.gz (as in BIDS) before sharing we can try to convert back to DICOM (research PACS suitable data). As much information is lost in the initial conversion this creates 'ugly' DICOM files.

```bash
./uglify -i data/ISLES-2022/sub-strokecase0001 -m data/ISLES-2022/derivatives/sub-strokecase0001 /tmp/bla
tree -L 3 /tmp/bla
/tmp/bla
└── strokecase0001
    └── 0001
        ├── 1_adc
        ├── 2_dwi
        ├── 3_FLAIR
        └── 4_msk
```

## Introduction

The ISLES-2022 dataset is an example volumetric medical imaging collection that contains images in the NIfTI format (extensions .nii.gz or .nii). Each .nii file is acompanied with a side-loading javascript object notation (json) text file that contains a single object with key-value pairs. Here an example:

```json
{
  "dataset": "ISLES_22",
  "ImageType": [
    "ORIGINAL",
    "PRIMARY",
    "M_IR",
    "M",
    "IR"
  ],
  "Modality": "MR",
  "Manufacturer": "Philips Medical Systems",
  "ManufacturerModelName": "Achieva dStream",
  "PatientSex": "F",
  "PatientAge": "088Y",
  "PatientWeight": 77.0,
  "BodyPartExamined": "BRAIN",
  "ScanningSequence": "IR",
  "SequenceVariant": "SK",
  "ScanOptions": "FS",
  "MRAcquisitionType": "3D",
  "SliceThickness": 1.0,
  "RepetitionTime": 4800.0,
  "EchoTime": 272.629,
  "InversionTime": 1650.0,
  "NumberOfAverages": 2.0,
  "ImagingFrequency": 127.76528,
  "ImagedNucleus": "1H",
  "EchoNumbers": 1,
  "MagneticFieldStrength": 3.0,
  "SpacingBetweenSlices": 0.712,
  "NumberOfPhaseEncodingSteps": 250,
  "EchoTrainLength": 170,
  "PercentSampling": 77.9164810180664,
  "PercentPhaseFieldOfView": 100.0,
  "PixelBandwidth": 1305.0,
  "ReconstructionDiameter": 250.0,
  "ReceiveCoilName": "MULTI COIL",
  "AcquisitionMatrix": [
    0,
    252,
    250,
    0
  ],
  "FlipAngle": 90.0,
  "PatientPosition": "HFS",
  "AcquisitionDuration": 235.20001220703125,
  "DiffusionBValue": 0.0,
  "DiffusionGradientOrientation": [
    0.0,
    0.0,
    0.0
  ],
  "ImagePositionPatient": [
    -35.706329450628,
    -112.80021829155,
    131.362080136887
  ],
  "ImageOrientationPatient": [
    -0.0277779083698,
    0.9990165233612,
    -0.0345597974956,
    0.05153298750519,
    -0.0330959893763,
    -0.9981227517127
  ],
  "Rows": 352,
  "Columns": 352
}
```

Most of these tags represent actual DICOM tags (not "dataset"). The goal of this project is to generate new DICOM files that are valid (for a given PACS system) and contain the image data (from the .nii) and the header information (from the json).

Warning: This approach is problematic because not all DICOM tags are present in the json file, only the once deemed useful by the original conversion tool (DICOM to .nii + .json) have been included. Our re-created DICOM files will therefore have missing information that might be required for data processing steps like minimally pre-processing. Even worse our generated DICOM files will pretend to be coming from a vendor without following their standard. A tool that tries to process or visualize such DICOM data might produce errors.

### Build

Use cmake and create either a 'Release' (fast) or a 'Debug' build.

```bash
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

To run a folder like ISLES:

```bash
for u in {0..250}; do 
   a=$(printf '%04d' $u); 
   ./uglify -i data/ISLES-2022/sub-strokecase${a} -m data/ISLES-2022/derivatives/sub-strokecase${a} /tmp/bla/
done
```