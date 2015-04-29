% % if running matlab 8.5:
%addpath(genpath('/a/software/matlab/8.5/toolbox/stats/'));

%run decision tree classification with 10-fold cross-validation

% load connectivity data features, 44 & 45 labels, and op & tri freesurfer labels
data = load('/scr/murg2/MachineLearning/BrocaDataAparc.mat');

% subselect data only for speed purposes
sel = [1:1000:size(data.features,2)]; % number of training sets
% sel_feat = [1:100:size(data.features,1)]; % number of conn features

features = [data.op(sel); data.tri(sel); data.features(:,sel)];
labels = data.labels(1,sel);
features = features';
labels = labels';

% run classification training
tree = ClassificationTree.fit(features, labels, 'crossval', 'on');
% For matlab 8.5 
% tree = fitctree(features, labels,'CrossVal','on');
% output cross-validation loss for each fold
L = kfoldLoss(tree,'mode','individual')
% visualize decision tree
view(tree.Trained{1}, 'Mode', 'graph');

%test classifier on a novel dataset
% import novel dataset
fid = fopen(['/scr/murg2/HCP500_glyphsets/994273/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');
% get freesurfer labels for new data
AnatLabels = gifti(['/scr/murg2/HCP500_glyphsets/994273/994273.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op_new = AnatLabelsData == 18;
tri_new = AnatLabelsData == 20;
op_new = op_new';
tri_new = tri_new';
data_new = [op_new; tri_new; M];
data_new = data_new';

% predict labels for novel dataset
label = predict(tree.Trained{1},data_new);
% save results for visualization
filename = ['/scr/murg2/MachineLearning/994273_results.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', label);
fclose(fid);