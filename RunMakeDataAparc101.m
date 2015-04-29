%runs MakeDataAparc script for 101 HCP subjects
%saves matrix with connectivity data (features), a vector of broca label
%data (labels), plus vectors of Op and Tri freesurfer labels
%NOTE: connectivity file is HUGE (15GB)! Don't try to open on a normal computer.

subject_list = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt']);

features = [];
labels = [];
op = [];
tri = [];

for i=1:length(subject_list)
[feat, lab, o, t] = MakeDataAparc(subject_list(i));

features = horzcat(features, feat);
labels = horzcat(labels, lab);
op = horzcat(op, o);
tri = horzcat(tri, t);

fprintf('Subject # %u\n',i)

end

clear feat
clear lab
clear o
clear t

save('BrocaDataAparc.mat', '-v7.3', 'features', 'labels', 'op', 'tri');