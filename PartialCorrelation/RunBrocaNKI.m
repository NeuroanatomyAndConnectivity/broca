subject_list = importdata(['/scr/murg2/MachineLearning/partialcorr/NKI/subject_list_NKI.txt']);

%% make group average time series and run ICA for NKI data

group_ts = zeros(10242,895);

for i=1:length(subject_list)
    ts = LoadNorm_ts_NKI(subject_list(i));
    group_ts = group_ts + ts;

    fprintf('Subject # %u\n',i)

end

group_ts = group_ts./length(subject_list);
group_ts_sparse = sparse(group_ts);

[ica10_ts, A_ts, W_ts] = fastica(group_ts_sparse, 'approach', 'symm', 'numOfIC', 10);
save('/scr/murg2/MachineLearning/partialcorr/NKI/ICA10_NKI.mat', 'ica10_ts', 'A_ts', 'W_ts');

%% create correlation matrices

for i=1:length(subject_list)
Write_CorrMat_NKI(subject_list(i));

fprintf('Subject # %u\n',i)

end

%% run Broca parcellation

for i=1:length(subject_list)
BrocaPartCorrICA10_NKI(subject_list(i));

fprintf('Subject # %u\n',i)

end

%% create group probability maps

prob45 = zeros(10242, 1);
prob44 = zeros(10242, 1);;

for i=1:length(subject_list)
    data = importdata(['/scr/murg2/MachineLearning/partialcorr/NKI/results/' char(subject_list(i)) '_partcorrMaxclust_rmtop2.1D']);
    ba44 = data;
    ba44(ba44==1)=0;
    ba44(ba44==2)=1;
    ba45 = data;
    ba45(ba45==2)=0;
    
    prob45 = prob45 + ba45;
    prob44 = prob44 + ba44;

    fprintf('Subject # %u\n',i)

end

prob44 = prob44./length(subject_list);
prob45 = prob45./length(subject_list);

filename = ['/scr/murg2/MachineLearning/partialcorr/NKI/results/GroupProb44_NKI.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob44);
fclose(fid);

filename = ['/scr/murg2/MachineLearning/partialcorr/NKI/results/GroupProb45_NKI.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob45);
fclose(fid);


