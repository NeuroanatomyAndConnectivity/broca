function[] = ClassifyBroca(subject)

subject = num2str(subject);

% get decision trees
tree45 = matfile('/scr/murg2/MachineLearning/DecisionTree_BA45.mat');
tree45 = tree45.tree_pruned;
tree44 = matfile('/scr/murg2/MachineLearning/DecisionTree_BA44.mat');
tree44 = tree44.tree_pruned;

% get connectivity from new dataset
fid = fopen(['/scr/murg2/HCP500_glyphsets/' subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% get pars opercularis and pars triangularis freesurfer labels for new dataset
AnatLabels = gifti(['/afs/cbs.mpg.de/documents/connectome/_all/' subject '/MNINonLinear/fsaverage_LR32k/' subject '.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
op = op';
tri = tri';
op = double(op);
tri = double(tri);

% get surface for new dataset
%may need to convert surface from .asc, use freesurfer mris_convert L.midthickness.asc lh.midthickness
surf = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.midthickness']);
surf = surfGetNeighbors(surf);

% calculate geodesic distance from pars opercularis
distOp = surfGeoDist_parcellation(surf, op);
% calculate geodesic distance from pars triangularis
distTri = surfGeoDist_parcellation(surf, tri);

% predict BA 45
mask = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast45-44.1D']);
mask = mask';
label45 = eval(tree45,M(:,find(mask)));
label45 = str2double(label45);
label45 = label45';
% apply anatomical constraint
results45 = label45;
results45(distTri>20)=0;
% save results
filename = ['/scr/murg2/MachineLearning/' subject '_results_BA45_constrained.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results45);
fclose(fid);

% predict BA 44
mask = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast44-45.1D']);
mask = mask';
label44 = eval(tree44,M(:,find(mask)));
label44 = str2double(label44);
label44 = label44';
% apply anatomical constraint
results44 = label44;
results44(distOp>10)=0;
% save results
filename = ['/scr/murg2/MachineLearning/' subject '_results_BA44_constrained.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results44);
fclose(fid);

%Combine BA 45 and 44 labels
results44(results44==1)=2;
results = results44 + results45;
results(results==3)=1.5;
filename = ['/scr/murg2/MachineLearning/' subject '_results_overlap_constrained.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);