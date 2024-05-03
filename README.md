# un-BIDS

For a lot of A.I. model training the publicly available data on kaggle etc. is useful. If the data has been converted to nii.gz (as in BIDS) before sharing we can try to convert back to DICOM (research PACS suitable data). As much information is lost in the initial conversion this creates 'ugly' DICOM files.

```bash
> ./uglify -i data/ISLES-2022/sub-strokecase0001 /tmp/bla
> tree -L 1 /tmp/bla
/tmp/bla
├── 1_adc
├── 2_dwi
└── 3_FLAIR
```

### Build

Use cmake and create either a 'Release' (fast) or a 'Debug' build.

```
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```
