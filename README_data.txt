Data created by RunMakeDataAparc101.m can be found on /scr/murg2/MachineLearning/BrocaDataAparc_new.mat
This file is huge (15GB)! It contains two variables:
features --> 32492x114718 matrix of connectivity data from 101 subjects within the vlpfc mask from the freesurfer tri + op labels
labels --> 1x114718 vector of corresponding labels where 1=BA45, 2=BA44, 0=neither

BrocaData_new.mat (5.5GB) is the same thing but without the zeros (so containing only data for manual BA 45 and 44 labels and no 'neither' regions)

List of subject IDs can be found at /scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt
