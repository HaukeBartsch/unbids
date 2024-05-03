# un-BIDS

For a lot of A.I. model training the publicly available data on kaggle etc. is useful. Sometimes those data have been converted from DICOM to nii.gz (as in BIDS). Here we try to undo this conversion to get data suitable for the research PACS. As much information is lost in the initial conversion this creates 'ugly' DICOM files.

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