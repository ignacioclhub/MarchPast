# MarchPastTesting

SETUP:
- Subjects: All subjects are from ADNI Multiband Database.
- Parcelation: Atlas used for parcelation is Schaefer2018_1000Parcels_17Networks. (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal).


FILES:
Scripts and single files needed:
- "test_MarchenkoPasture.m": general script which plots several testings. Feel free to improve or change it.
- "preparemarchenko.m": function to compute all what MP needs.
- "Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.mat": File with parcelation labels and coordinates.
- "Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii": MRI atlas with Parcelation labels.
Data:
- "/Data/Concatenated/" - Folder containing 3 files, 1 x group (AD-Alzheimer disease,CN-Healthy Subjects, MCI-Mild cognitive impairment). Each file contains 10 subjects concatenated timeseries (9760timex1000parcels).
- "/Data/Raw/" - Folder containing 3 folders, 1 x group (AD-Alzheimer disease,CN-Healthy Subjects, MCI-Mild cognitive impairment). Each folder contains 10 subjects timeseries (976timex1000parcels).
