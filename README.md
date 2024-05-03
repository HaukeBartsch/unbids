# un-BIDS

For a lot of A.I. model training the publicly available data on kaggle etc. is useful. Sometimes those data have been converted from DICOM to nii.gz (as in BIDS). Here we try to undo this conversion to get data suitable for the research PACS. As much information is lost in the initial conversion this creates 'ugly' DICOM files.



### Build

Use cmake and create either a 'Release' (fast) or a 'Debug' build.

```
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```