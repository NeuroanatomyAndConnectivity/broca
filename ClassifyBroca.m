function[cost] = ClassifyBroca(subject)

subject = num2str(subject);

%% Classify new dataset using decision tree

% get decision trees
tree = matfile('/scr/murg2/MachineLearning/DecisionTree4445_AvgConnFeat_unpruned.mat');
tree = tree.tree;

% get new dataset
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

% calculate geodesic distances
distOp = surfGeoDist_parcellation(surf, op);
distTri = surfGeoDist_parcellation(surf, tri);

mask44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast44-45.1D']);
mask45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast45-44.1D']);
mask = mask44 + mask45;
mask = mask';
mask(isnan(mask))=0;
sel_feat = find(mask);
data_new = [distOp; distTri; M(sel_feat,:)];
data_new = data_new';

label = eval(tree, data_new);
cost = test(tree, 'test', data_new, label);

label = str2double(label);
label = label';

%% apply spatial constraint according to distance from op and tri

optri = op + tri;
distOpTri = surfGeoDist_parcellation(surf, optri);
label(distOpTri>35)=0;
%% apply spatial constraint using pars orbitalis as a hard cut-off

orb = AnatLabelsData == 19;
orb = orb';
orb = double(orb);
label(orb==1)=0;

%% apply spatial constraint according to neighbors

surf.nbr(surf.nbr==0)=1;
results = label;
for i = 1:length(surf.nbr)
    col = surf.nbr(:,i);
    col2 = results(col);
    sumcol0 = sum(col2(:) == 0);
    sumcol1 = sum(col2(:) == 1);
    sumcol2 = sum(col2(:) == 2);
    if sumcol1>3
    results(:,i)=1;
    elseif sumcol2>3
    results(:,i)=2;
    end
    if sumcol0>4
    results(:,i)=0;
    end
end

%% apply spatial constraint using cluster size

% create edge list for BA45
label45 = results;
label45(label45==2) = 0;
edgelist45 = [];
for i = 1:length(label45)
    nbrs = surf.nbr(:,i);
    nbrslabel = label45(nbrs);
    for j = 1:length(nbrs)
        if label45(i)~=0 && label45(i) == nbrslabel(j)
            edge = [i nbrs(j)];
            edgelist45 = vertcat(edgelist45, edge);
        else
            edgelist45 = edgelist45;
        end
    end
end

% build adjacency matrix from edge list for BA45
A45 = sparse([edgelist45(:,1),edgelist45(:,2)],[edgelist45(:,2),edgelist45(:,1)],1);

% get network components for BA45
[nComponents45,sizes45,members45] = networkComponents(A45);
% find the largest component of BA45
largest45 = members45{1};
results45 = zeros(size(results));
results45(largest45) = 1;

% create edge list for BA44
label44 = results;
label44(label44==1) = 0;
edgelist44 = [];
for i = 1:length(label44)
    nbrs = surf.nbr(:,i);
    nbrslabel = label44(nbrs);
    for j = 1:length(nbrs)
        if label44(i)~=0 && label44(i) == nbrslabel(j)
            edge = [i nbrs(j)];
            edgelist44 = vertcat(edgelist44, edge);
        else
            edgelist44 = edgelist44;
        end
    end
end

% build adjacency matrix from edge list for BA44
A44 = sparse([edgelist44(:,1),edgelist44(:,2)],[edgelist44(:,2),edgelist44(:,1)],1);

% get network components for BA44
[nComponents44,sizes44,members44] = networkComponents(A44);
% find the largest component of BA44
largest44 = members44{1};
results44 = zeros(size(results));
results44(largest44) = 2;

%% combine results from BA 44 and 45 and save

results = results45 + results44;

filename = ['/scr/murg2/MachineLearning/newdata_results/' subject '_unpruned_constraint_clust.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);