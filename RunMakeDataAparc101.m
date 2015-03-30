subject_list = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt']);

Data = [];

for i=1:length(subject_list)
data = MakeDataAparc(subject_list(i));

Data = horzcat(Data, data);

fprintf('Subject # %u\n',i)

end

save('BrocaDataAparc.m', '-v7.3', 'Data');