function[corrmat] = Write_CorrMat_NKI(subject)

subject = char(subject);

tseries = h5read(['/scr/ilz1/nilearn_vol2surf_sink/fwhm6/' subject '_lh_preprocessed_fsaverage5_fwhm6.hdf5'],'/mat');
corrmat = corrcoef(tseries);

fileID = fopen(['/scr/murg2/NKI/' subject '/' subject '_lh_preprocessed_fsaverage5_fwhm6_corr.bin'],'w');
fwrite(fileID, corrmat, 'single');
fclose(fileID);