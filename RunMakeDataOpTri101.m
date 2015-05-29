%runs MakeDataOpTri script for 101 HCP subjects

subject_list = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt']);

op = [];
tri = [];

for i=1:length(subject_list)
[Op, Tri] = MakeDataOpTri(subject_list(i));

op = horzcat(op, Op');
tri = horzcat(tri, Tri');

fprintf('Subject # %u\n',i)

end

save('BrocaDataOpTri.mat', '-v7.3', 'op', 'tri');