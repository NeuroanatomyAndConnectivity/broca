%runs MakeDataAparc script for 101 HCP subjects
%saves matrix with connectivity data (features) and a vector of broca label data (labels) separately
%NOTE: connectivity file is HUGE (15GB)! Don't try to open on a normal computer.

subject_list = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt']);

Data = [];

for i=1:length(subject_list)
data = MakeDataAparc(subject_list(i));

Data = horzcat(Data, data);

fprintf('Subject # %u\n',i)

end

labels = Data(32493,:);
features = Data(1:32492,:);

save('BrocaDataAparc.mat', '-v7.3', 'labels', 'features');