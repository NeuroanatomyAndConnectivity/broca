function[] = BrocaPartCorr(subject)

subject = num2str(subject);

%% Import data and run partial correlation for each area

% get correlation matrix for individual subject
fid = fopen(['/scr/murg1/HCP500_glyphsets/' subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% get group connectivity maps
conn44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast44-45.1D']);
conn44(isnan(conn44))=0;
conn45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_thres_contrast45-44.1D']);
conn45(isnan(conn45))=0;

% get pars opercularis and pars triangularis freesurfer labels for new dataset
AnatLabels = gifti(['/afs/cbs.mpg.de/documents/connectome/_all/' subject '/MNINonLinear/fsaverage_LR32k/' subject '.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
op = op';
tri = tri';
op = double(op);
tri = double(tri);

% import manual probability maps and binarize them
prob44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_mean_101.1D']);
prob45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_mean_101.1D']);
prob44(prob44>0)=1;
prob45(prob45>0)=1;

% add freesurfer labels to manual prob maps to create broca mask
broca = prob44' + prob45' + op + tri;

% mask the correlation matrix to Broca
M_small = M(:,find(broca));
M_small(isnan(M_small))=0;

% run partial correlation for 45 vs 44
partcorr45 = partialcorr(M_small, conn45, conn44);

% run partial correlation for 44 vs 45
partcorr44 = partialcorr(M_small, conn44, conn45);

% make results wholebrain and save as 1D
results45 = zeros(length(M),1);
results45(find(broca)) = partcorr45;
filename = ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorr45.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results45);
fclose(fid);

results44 = zeros(length(M),1);
results44(find(broca)) = partcorr44;
filename = ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorr44.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results44);
fclose(fid);

%% vizualize results
% may need to convert surface from .asc, use freesurfer mris_convert L.midthickness.asc lh.midthickness
surf = SurfStatReadSurf1(['/scr/murg1/HCP500_glyphsets/' subject '/lh.very_inflated']);
figure('visible','off'); SurfStatView(results45,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorr45.png']); close all;
figure('visible','off'); SurfStatView(results44,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorr44.png']); close all;

%% Apply spatial constraint using geodesic distance

% calculate geodesic distances from freesurfer labels
surf = SurfStatReadSurf1(['/scr/murg1/HCP500_glyphsets/' subject '/lh.midthickness']);
surf = surfGetNeighbors(surf);
distOp = surfGeoDist_parcellation(surf, op);
distTri = surfGeoDist_parcellation(surf, tri);

%invert and then normalize geodist maps
invOp = distOp*-1;
invmin = min(invOp);
invmax = max(invOp);
invrange = invmax - invmin;
normdistOp = (invOp - invmin) / invrange;

invTri = distTri*-1;
invmin = min(invTri);
invmax = max(invTri);
invrange = invmax - invmin;
normdistTri = (invTri - invmin) / invrange;

% multiply partial correlation results by inverse distance maps
new45 = results45.*normdistTri';
new44 = results44.*normdistOp';

% normalize partcorr maps
min45 = min(new45);
max45 = max(new45);
range45 = max45 - min45;
norm45 = (new45 - min45) / range45;
min44 = min(new44);
max44 = max(new44);
range44 = max44 - min44;
norm44 = (new44 - min44) / range44;

% % create "neither" map by adding 44 and 45 correlations and inverting
% neither = (new45 + new44)*-1;
% or import neither prob map
neither = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/neitherWholebrain_mean.1D']);

% apply winner-take-all to create binary partition from the partial correlation maps
maps = [norm45, norm44, neither];
[val,part] = max(maps,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% save results
surf = SurfStatReadSurf1(['/scr/murg1/HCP500_glyphsets/' subject '/lh.very_inflated']);
figure('visible','off'); SurfStatView(part,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorrResults.png']); close all;

part(part==3) = 0;
filename = ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorrResults.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', part);
fclose(fid);


%% apply spatial constraint using cluster size

% get neighborhood information
surf = SurfStatReadSurf1(['/scr/murg1/HCP500_glyphsets/' subject '/lh.very_inflated']);
surf = surfGetNeighbors(surf);

surf.nbr(surf.nbr==0)=1; %necessary for correct indexing

% create edge list for BA45
label45 = part';
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
results45 = zeros(size(part));
results45(largest45) = 1;

% create edge list for BA44
label44 = part;
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
results44 = zeros(size(part));
results44(largest44) = 2;

% combine results from BA 44 and 45
results = results45 + results44;
filename = ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorrResults_maxclust.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);

figure('visible','off'); SurfStatView(results,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorrResults_maxclust.png']); close all;