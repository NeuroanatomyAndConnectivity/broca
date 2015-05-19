addpath(genpath('/scr/murg2/scripts/'));
cd /scr/murg2/scripts/FEAST/
CompileFEAST

% load connectivity data features, 44 & 45 labels, and op & tri freesurfer labels
fprintf('loading data...')
data = load('/scr/murg2/MachineLearning/BrocaDataAparc.mat');

% subselect data only for speed purposes
fprintf('subselecting data...')
sel = [1:2:size(data.features,2)]; % number of training sets
%sel_feat = [1:10:size(data.features,1)]; % number of conn features

%features = [data.op(sel); data.tri(sel); data.features(:,sel)];
features = data.features(:,sel);
labels = data.labels(1,sel);
features = features';
labels = labels';
op = data.op;
tri = data.tri;

fprintf('clearing some memory...')
clear data

%feature selection
fprintf('selecting features...')
features(isnan(features))=0;
selected = feast('jmi', 500, features, labels); %adjust number of features to select

%save selected nodes as 1D
nodes = zeros(1,32492);
nodes(1,selected)=1;
filename = ['/scr/murg2/MachineLearning/selected_nodes_everysecond.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', nodes);
fclose(fid);
    
% run classification training
fprintf('running classification...')
tree = classregtree(features(:,selected), labels, 'method', 'classification', 'minparent', 10)

% cross-validate decision tree with 10-fold
[cost,secost,ntnodes,bestlevel] = test(tree,'crossvalidate', features(:,selected), labels)
%prune decision tree and then view
tree_pruned = prune(tree,'level',bestlevel);
view(tree_pruned)

%save pruned decision tree
save('DecisionTree500Features_everysecond.mat', '-v7.3', 'tree_pruned', 'selected');

%test classifier on a novel dataset

% import novel dataset
fid = fopen(['/scr/murg2/HCP500_glyphsets/984472/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% predict labels for novel dataset
label = eval(tree_pruned,M(:,selected));
label = str2double(label);
% save results for visualization
filename = ['/scr/murg2/MachineLearning/984472_results_500feat.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', label);
fclose(fid);

% % or if freesurfer labels were included in training:
% 
% % get freesurfer labels for new data
% AnatLabels = gifti(['/scr/murg2/HCP500_glyphsets/984472/984472.L.aparc.32k_fs_LR.label.gii']);
% AnatLabelsData = AnatLabels.cdata;
% op_new = AnatLabelsData == 18;
% tri_new = AnatLabelsData == 20;
% op_new = op_new';
% tri_new = tri_new';
% data_new = [op_new; tri_new; M];
% data_new = data_new';
% 
% % predict labels for novel dataset
% label = eval(tree,data_new);
% label = str2double(label);
% % save results for visualization
% filename = ['/scr/murg2/MachineLearning/984472_results_new_5feat.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', label);
% fclose(fid);