%% get subject list

subject_list = importdata(['/scr/murg2/MachineLearning/partialcorr/NKI/subject_list_NKI.txt']);

%% create correlation matrices

for i=1:length(subject_list)
Write_CorrMat_NKI(subject_list(i));

fprintf('Subject # %u\n',subject_list(i))

end

%% run Broca parcellation

for i=1:length(subject_list)
    try
    BrocaPartCorrNKI(subject_list(i), 0.1, 'rm_0p1');
%     BrocaPartCorrNKI(subject_list(i), 0.2, 'rm_0p2');
%     BrocaPartCorrNKI(subject_list(i), 0.3, 'rm_0p3');
%     BrocaPartCorrNKI(subject_list(i), 0.4, 'rm_0p4');
%     BrocaPartCorrNKI(subject_list(i), 0.5, 'rm_0p5');
    catch
        fprintf('Could not run subject # %u\n', i)
        continue
    end
    fprintf('Subject # %u\n', i)
end


%% create group probability maps

prob45 = zeros(10242, 1);
prob44 = zeros(10242, 1);
rm = [];

for i=1:length(subject_list)
    try
    data = importdata(['/scr/murg2/MachineLearning/partialcorr/20comps_results/NKI_ICA_indiv/' char(subject_list(i)) '_ICA_indiv_SW_rm_0p1.1D']);
    ba44 = data;
    ba44(ba44==1)=0;
    ba44(ba44==2)=1;
    ba45 = data;
    ba45(ba45==2)=0;
    
    prob45 = prob45 + ba45;
    prob44 = prob44 + ba44;
    catch
        rm = horzcat(rm, i);
        fprintf('Could not import subject # %u\n',i)
        continue
    end
    fprintf('Subject # %u\n', i)

end

subject_list(rm, :) = [];

prob44 = prob44./length(subject_list);
prob45 = prob45./length(subject_list);

filename = ['/scr/murg2/MachineLearning/partialcorr/20comps_results/NKI_ICA_indiv/GroupProb44_NKI_ICA_indiv_SW_rm_0p1.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob44);
fclose(fid);

filename = ['/scr/murg2/MachineLearning/partialcorr/20comps_results/NKI_ICA_indiv/GroupProb45_NKI_ICA_indiv_SW_rm_0p1.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob45);
fclose(fid);

%% evaluate spatial overlap for 8 manually labeled datasets
subject_list = subject_list(1:8);
for i=1:length(subject_list)
    try
    [~, ~, DC(i,1)] = AutoVsManualDC_NKI(subject_list(i), '_ICA_indiv_SW_rm_0p1');
    [~, ~, DC(i,2)] = AutoVsManualDC_NKI(subject_list(i), '_ICA_indiv_SW_rm_0p2');
    [~, ~, DC(i,3)] = AutoVsManualDC_NKI(subject_list(i), '_ICA_indiv_SW_rm_0p3');
    [~, ~, DC(i,4)] = AutoVsManualDC_NKI(subject_list(i), '_ICA_indiv_SW_rm_0p4');
    [~, ~, DC(i,5)] = AutoVsManualDC_NKI(subject_list(i), '_ICA_indiv_SW_rm_0p5');
    catch
        fprintf('Could not run subject # %u\n',i)
        continue
    end
    fprintf('Subject # %u\n', i)

end
