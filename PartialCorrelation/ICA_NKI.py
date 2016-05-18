import nibabel as nib
import numpy as np
from scipy.io import savemat
from nilearn.decomposition.canica import CanICA
from nibabel import gifti

# LEFT indices
img = nib.load('/a/documents/connectome/_all/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii')
cort = img.header.matrix.mims[1].brainModels[0].vertexIndices.indices

# load data
#subs = ['A00028185'] # list of subject ids to include
subs = np.loadtxt("/scr/murg2/MachineLearning/partialcorr/NKI/subject_list_NKI.txt", dtype=str)
filenames = []
for sub in subs:
    print "SUBJECT ", sub
    try:
        data = gifti.giftiio.read('/scr/ilz1/nilearn_vol2surf_sink/fwhm6/%s_lh_preprocessed_fsaverage5_fwhm6.gii' % sub)
        data = data[cort, :]
        data = np.expand_dims(data,axis=1)    
        data = np.expand_dims(data,axis=1)
        filename = '/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/%s.nii.gz' % sub
        img = nib.Nifti1Image(data, np.eye(4))
        img.to_filename(filename)
        filenames.append(filename)
    except:
        print "subject " + sub + " cannot be run"
         
    # load data:                
    # data = nib.load('/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/HCP_Q1-Q6_R468_rfMRI_groupPCA_d4500_Eigenmaps_MGTR_left_metric_smoothed2.dtseries.nii')
    # d = data.get_data().squeeze()
    # cort = np.where(d.sum(axis=0) != 0.)[0]
    # d = d[:,cort]
    # d = np.expand_dims(d.T,axis=1)
    #d = np.expand_dims(d,axis=1)
    #filename = '%s.nii.gz' % sub
    #img = nib.Nifti1Image(d, np.eye(4))
    #img.to_filename(filename)

# create artificial mask:
mask = np.ones(29696)
mask = np.expand_dims(mask,axis=1)
mask = np.expand_dims(mask,axis=1)
img = nib.Nifti1Image(mask, np.eye(4))
img.to_filename('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/mask.nii.gz')

# run ICA on group level:
n_components = 10
n_jobs = 10 #number of CPUs used
canica = CanICA(mask='/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/mask.nii.gz',
                n_components=n_components, smoothing_fwhm=0.,
                threshold=None, verbose=10, random_state=0, n_jobs=n_jobs)

canica.fit(filenames)

# Retrieve the independent components in brain space
components_img = canica.masker_.inverse_transform(canica.components_)

A = np.zeros((32492,n_components))
# STILL NEED TO RESOLVE CORT INDICIES

# LEFT indices
img = nib.load('/a/documents/connectome/_all/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii')
cort = img.header.matrix.mims[1].brainModels[0].vertexIndices.indices

A[cort,:] = components_img.get_data().squeeze()
#np.save('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_HCP/ica_HCP500_output_%s.npy' % str(n_components), A)
savemat('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/ica_NKI_output_%s.mat' % str(n_components), {'ic':A})

## RUN ICA ON INDIVIDUAL LEVEL:
for sub in subs:
    try:
        n_components = 10
        n_jobs = 10 #number of CPUs used
        canica = CanICA(mask='/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/mask.nii.gz',
                        n_components=n_components, smoothing_fwhm=0.,
                        threshold=None, verbose=10, random_state=0, n_jobs=n_jobs)
        filename = '/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/%s.nii.gz' % sub
        canica.fit(filename)
        
        # Retrieve the independent components in brain space
        components_img = canica.masker_.inverse_transform(canica.components_)
        
        A = np.zeros((32492,n_components))
        # STILL NEED TO RESOLVE CORT INDICIES
        A[cort,:] = components_img.get_data().squeeze()
        #np.save('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/ica_output_%s_%s.npy' % (sub, str(n_components)), A)
        savemat('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_NKI/ica_output_%s_%s.mat' % (sub, str(n_components)), {'ic':A})
    except:
        print "subject " + sub + " cannot be run"
        
        
        
        
