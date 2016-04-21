subject_list = importdata(['/scr/murg2/MachineLearning/partialcorr/NKI/subject_list_NKI.txt']);

for i=1:length(subject_list)
Write_CorrMat_NKI(subject_list(i));

fprintf('Subject # %u\n',i)

end

for i=1:length(subject_list)
BrocaPartCorrICA10_NKI(subject_list(i));

fprintf('Subject # %u\n',i)

end