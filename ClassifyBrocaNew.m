function[cost] = ClassifyBrocaNew(subject)

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
% may need to convert surface from .asc, use freesurfer mris_convert L.midthickness.asc lh.midthickness
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

% import manual probability maps and binarize them
prob44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_mean_101.1D']);
prob45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_mean_101.1D']);
prob44(prob44>0)=1;
prob45(prob45>0)=1;
% add freesurfer labels
broca = prob44' + prob45' + op + tri;

% mask the input data so only broca gets classified
data_new = data_new(:,find(broca));
data_new = data_new';

% run classification
[label, nodes] = eval(tree, data_new);
cost = test(tree, 'test', data_new, label);
label = str2double(label);
label = label';

% make results same size as surface
results = zeros([length(broca) 1]);
results(find(broca)) = label;

% filename = ['/scr/murg2/MachineLearning/testing/' subject '_unpruned_masked.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', results);
% fclose(fid);

% %view results in surfstat
% surf2 = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
% SurfStatViewData(results, surf2)

% make probability maps for each area and then save
probmap44 = zeros(size(label));
probmap45 = zeros(size(label));
probmapNeither = zeros(size(label));

prob_class = classprob(tree);
for i = 1:length(label)
        nodenum = nodes(i);
        probmap45(i)=prob_class(nodenum,2);
        probmap44(i)=prob_class(nodenum,3);
        probmapNeither(i)=prob_class(nodenum,1);
end

BA44prob = zeros([length(broca) 1]);
BA44prob(find(broca)) = probmap44;

% filename = ['/scr/murg2/MachineLearning/testing/' subject '_class_probmap44.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', BA44prob);
% fclose(fid);
% 
BA45prob = zeros([length(broca) 1]);
BA45prob(find(broca)) = probmap45;

% filename = ['/scr/murg2/MachineLearning/testing/' subject '_class_probmap45.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', BA45prob);
% fclose(fid);
% 
Neitherprob = zeros([length(broca) 1]);
Neitherprob(find(broca)) = probmapNeither;

% filename = ['/scr/murg2/MachineLearning/testing/' subject '_class_probmapNeither.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', Neitherprob);
% fclose(fid);
% 
% % %% apply spatial constraint according to distance from op and tri
% % 
% optri = op + tri;
% distOpTri = surfGeoDist_parcellation(surf, optri);
% label(distOpTri>35)=0;

%% apply spatial constraint according to probability of neighbors

results = results';
surf.nbr(surf.nbr==0)=1;
for i = 1:length(surf.nbr)
    nbrs = surf.nbr(:,i);
    nbrs44 = BA44prob(nbrs);
    nbrs45 = BA45prob(nbrs);
    nbrsNeither = Neitherprob(nbrs);
    sum44 = nansum(nbrs44);
    sum45 = nansum(nbrs45);
    sumNeither = nansum(nbrsNeither);
    sum = [sumNeither; sum45; sum44];
    [M,I] = max(sum);
    results(:,i)=I-1;
end

% filename = ['/scr/murg2/MachineLearning/testing/' subject '_class_prob_nbrs.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', results);
% fclose(fid);

%% apply spatial constraint using pars orbitalis as a hard cut-off

orb = AnatLabelsData == 19;
orb = orb';
orb = double(orb);
results(orb==1)=0;

%% clean up results using morphological operations on adjacency matrix

% surf.nbr(surf.nbr==0)=1; %necessary for correct indexing
% results = results';

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
A45_full = full(A45);

%define a structuring element for the neighborhood
NH = strel('disk', 10) %disk-shaped neighborhood with radius of N

%perform morphological operations using strel neighborhood definition
A45_clean = imclose(A45_full,NH);

%clean up the results using bwmorph
%NB: turns adjacency matrix from double to logical
A45_clean = bwmorph(A45_clean, 'bridge', Inf); %bridge gaps
% A45_clean = bwmorph(A45_clean, 'fill', Inf); %fill in isolated holes
A45_clean = bwmorph(A45_clean, 'clean', Inf); %remove isolated points
A45_clean = imfill(A45_clean,'holes');

% % turn the adjacency matrix back into a label vector
% v = sum(A45_clean);
% index_list = find(v ~= 0);
% results45 = zeros(size(results));
% results45(index_list) = 1;

%% take only the largest component

% find the largest component of BA45
[nComponents45,sizes45,members45] = networkComponents(A45_clean);
largest45 = members45{1};
results45 = zeros(size(results));
results45(largest45) = 1;

% %view results
% surf2 = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
% figure;
% SurfStatViewData(label45, surf2)
% figure;
% SurfStatViewData(results45, surf2)


%% repeat for BA44...
%% clean up results using morphological operations on adjacency matrix

% surf.nbr(surf.nbr==0)=1; %necessary for correct indexing
% results = results';

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
A44_full = full(A44);

%define a structuring element for the neighborhood
NH = strel('disk', 10) %disk-shaped neighborhood with radius of N

%perform morphological operations using strel neighborhood definition
A44_clean = imclose(A44_full,NH);

%clean up the results using bwmorph
%NB: turns adjacency matrix from double to logical
A44_clean = bwmorph(A44_clean, 'bridge', Inf); %bridge gaps
% A44_clean = bwmorph(A44_clean, 'fill', Inf); %fill in isolated holes
A44_clean = bwmorph(A44_clean, 'clean', Inf); %remove isolated points
A44_clean = imfill(A44_clean,'holes');

% % turn the adjacency matrix back into a label vector
% v = sum(A44_clean);
% index_list = find(v ~= 0);
% results44 = zeros(size(results));
% results44(index_list) = 1;

%% take only the largest component

% find the largest component of BA44
[nComponents44,sizes44,members44] = networkComponents(A44_clean);
largest44 = members44{1};
results44 = zeros(size(results));
results44(largest44) = 2;

% %view results
% surf2 = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
% figure;
% SurfStatViewData(label44, surf2)
% figure;
% SurfStatViewData(results44, surf2)

%% combine results from BA 44 and 45 and save

results_clean = results45 + results44;

surf2 = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
figure;
SurfStatViewData(results, surf2)
figure;
SurfStatViewData(results_clean, surf2)

filename = ['/scr/murg2/MachineLearning/testing/' subject '_class_prob_clean.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results_clean);
fclose(fid);