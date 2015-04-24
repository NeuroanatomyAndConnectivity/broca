%run decision tree classification with 10-fold cross-validation
load /scr/murg2/MachineLearning/BrocaData_new.mat
tree = ClassificationTree.fit(features', labels', 'crossval', 'on')
view(tree.Trained{1}, 'Mode', 'graph')

%import novel dataset and mask to freesurfer vlpfc label
AnatLabels = gifti(['/afs/cbs.mpg.de/documents/connectome/_all/984472/MNINonLinear/fsaverage_LR32k/984472.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
vlpfc = op + tri;
fid = fopen(['/scr/murg2/HCP500_glyphsets/984472/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');
M_masked = M(find(vlpfc),:);

%predict labels for novel dataset
label = predict(tree.Trained{1},M_masked);

%save results for visualization
results = zeros([length(vlpfc) 1]);
results(find(vlpfc)) = label;
filename = ['984472_results.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);