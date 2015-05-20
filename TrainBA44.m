addpath(genpath('/scr/murg2/scripts/'));

% load connectivity data features, 44 & 45 labels, and op & tri freesurfer labels
fprintf('loading data...')
data = load('/scr/murg2/MachineLearning/BrocaDataBA44.mat');

features = data.BA44_features;
labels = data.BA44_labels;
features = features';
labels = labels';
op = data.Op;

% only include features from contrasted group connectivity maps
mask = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast44-45.1D']);
mask = mask';
features = features(:,find(mask));
features(isnan(features))=0;

fprintf('clearing some memory...')
clear data

% run classification training
fprintf('running classification training...')
tree = classregtree(features, labels, 'method', 'classification', 'minparent', 10)

% cross-validate decision tree with 10-fold
fprintf('cross-validating...')
[cost,secost,ntnodes,bestlevel] = test(tree,'crossvalidate', features, labels);
%prune decision tree and then view
fprintf('pruning decision tree...')
tree_pruned = prune(tree,'level',bestlevel);
view(tree_pruned)

%save pruned decision tree
fprintf('saving decision tree...')
save('DecisionTree_BA44.mat', '-v7.3', 'tree_pruned');