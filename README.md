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