%% make group average time series and run ICA for NKI data

subject_list = importdata(['/scr/murg2/MachineLearning/partialcorr/NKI/subject_list_NKI.txt']);

group_ts = zeros(10242,895);

for i=1:length(subject_list)
    ts = LoadNorm_ts_NKI(subject_list(i));
    group_ts = group_ts + ts;

    fprintf('Subject # %u\n',i)

end

group_ts = group_ts./length(subject_list);
% figure; plot(group_ts(1,:));

group_ts_sparse = sparse(group_ts);

[ica10_ts, A_ts, W_ts] = fastica(group_ts_sparse, 'approach', 'symm', 'numOfIC', 10);
save('/scr/murg2/MachineLearning/partialcorr/NKI/ICA10_NKI.mat', 'ica10_ts', 'A_ts', 'W_ts');

surf = SurfStatReadSurf1('/scr/murg2/MachineLearning/partialcorr/NKI/lh.inflated');
figure; SurfStatView((A_ts(:,6)),surf);
