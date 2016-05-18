%% get subject list

subject_list = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/subject_list_101.txt']);

% subject_list = importdata(['/scr/murg2/MachineLearning/EightNewSubjects.txt']);
% subject_list = importdata(['/scr/murg1/HCP500_glyphsets/subject_list_HCP500.txt']);

%% run Broca parcellation

for i=1:length(subject_list)
    try
    BrocaPartCorrICAIndiv(subject_list(i), 0.4, 'rm_abs0p4');
    catch
        fprintf('Could not run subject # %u\n',subject_list(i))
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))
end

%% evaluate spatial overlap of manual and automated labels

for i=1:length(subject_list)
    try
    [DC(i,1), DC(i,2), DC(i,3)] = AutoVsManualDC(subject_list(i), '_ICA_indiv_SW_rm_abs0p4');
    catch
        fprintf('Could not run subject # %u\n',i)
        continue
    end
    fprintf('Subject # %u\n',subject_list(i))

end

%% create group probability maps

prob45 = zeros(32492, 1);
prob44 = zeros(32492, 1);
rm = [];

for i=1:length(subject_list)
    try
    data = importdata(['/scr/murg2/MachineLearning/partialcorr/20comps_results/HCP_ICA_Indiv/' num2str(subject_list(i)) '_ICA_indiv_SW_rm_abs0p4.1D']);
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
    fprintf('Subject # %u\n',subject_list(i))

end

subject_list(rm, :) = [];

prob44 = prob44./length(subject_list);
prob45 = prob45./length(subject_list);

filename = ['/scr/murg2/MachineLearning/partialcorr/20comps_results/HCP_ICA_Indiv/GroupProb44_HCP101_ICA_indiv_SW_rm_abs0p4.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob44);
fclose(fid);

filename = ['/scr/murg2/MachineLearning/partialcorr/20comps_results/HCP_ICA_Indiv/GroupProb45_HCP101_ICA_indiv_SW_rm_abs0p4.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', prob45);
fclose(fid);