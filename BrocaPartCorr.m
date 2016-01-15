function[] = BrocaPartCorr(subject)

subject = num2str(subject);

%% Import data and run partial correlation for each area

% get correlation matrix for individual subject
fid = fopen(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% get group connectivity maps
conn44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_nothresh_44.1D']);
conn44(isnan(conn44))=0;
conn45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_nothresh_45.1D']);
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
broca(broca>0)=1;

% mask the correlation matrix to Broca
M_small = M(:,find(broca));
M_small(isnan(M_small))=0;

% import ICA maps
ica = load('/scr/murg2/MachineLearning/partialcorr/ICA/ica_output_20.mat');
ica = ica.ic;
surf = SurfStatReadSurf1(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/lh.inflated']);

% remove components that look most like 45
ica_rm45 = ica;
ica_rm45(:,[1, 4, 10, 12, 17, 19]) = [];

% remove components that look most like 44
ica_rm44 = ica;
ica_rm44(:,[5, 9, 11, 15]) = [];

% run partial correlation for 45
partcorr45 = partialcorr(M_small, conn45, ica_rm45);

% run partial correlation for 44
partcorr44 = partialcorr(M_small, conn44, ica_rm44);

%run partial correlation for remaining ICA components
ica_rmboth = ica;
ica_rmboth(:,[1, 4, 5, 9, 10, 11, 12, 15, 17, 19]) = [];
partcorrIC = [];
for i=1:size(ica_rmboth,2)
    rm = ica_rmboth;
    rm(:,i) = [];
    rm = horzcat(conn45, conn44, rm);
    partcorr = partialcorr(M_small, ica_rmboth(:,i), rm);
    results = zeros(length(M),1);
    results(find(broca)) = partcorr;
    partcorrIC = horzcat(partcorrIC, results);
end

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

%% apply spatial contraint using geodesic distance calculated from...

% % % freesurfer labels
% % surf = SurfStatReadSurf1(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/lh.midthickness']);
% % surf = surfGetNeighbors(surf);
% % distOp = surfGeoDist_parcellation(surf, op);
% % distTri = surfGeoDist_parcellation(surf, tri);
% % invOp = distOp*-1;
% % invmin = min(invOp);
% % invmax = max(invOp);
% % invrange = invmax - invmin;
% % normdistOp = (invOp - invmin) / invrange;
% % invTri = distTri*-1;
% % invmin = min(invTri);
% % invmax = max(invTri);
% % invrange = invmax - invmin;
% % normdistTri = (invTri - invmin) / invrange;
% % new45 = results45.*normdistTri';
% % new44 = results44.*normdistOp';
% 
% % OR either max of partial correlation maps or max of manual probability maps
% % max44 = results44;
% max44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_mean_101.1D']);
% % max45 = results45;
% max45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_mean_101.1D']);
% max44(max44<(max(max44)))=0;
% max44(max44>0)=1;
% max45(max45<(max(max45)))=0;
% max45(max45>0)=1;
% 
% surf = SurfStatReadSurf1(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/lh.midthickness']);
% surf = surfGetNeighbors(surf);
% distmax44 = surfGeoDist_parcellation(surf, max44');
% distmax45 = surfGeoDist_parcellation(surf, max45');
% 
% % mask, invert, and then normalize geodist maps
% inv44 = distmax44.*broca;
% inv44 = inv44*-1;
% normdist44 = (inv44 - min(inv44)) / (max(inv44) - min(inv44));
% % normdist44 = normdist44*10;
% 
% inv45 = distmax45.*broca;
% inv45 = inv45*-1;
% normdist45 = (inv45 - min(inv45)) / (max(inv45) - min(inv45));
% % normdist45 = normdist45*10;
% 
% % multiply partial correlation results by inverse distance maps
% new45 = results45.*normdist45';
% new44 = results44.*normdist44';

% % OR use geodist from thresholded broca prob map 
% manprob = importdata(['/scr/murg2/MachineLearning/partialcorr/GroupProbBroca.1D']);
% manprob(manprob>=0.5)=1;
% manprob(manprob<1)=0;
% surf = SurfStatReadSurf1(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/lh.midthickness']);
% surf = surfGetNeighbors(surf);
% distbroca = surfGeoDist_parcellation(surf, manprob');
% % mask, invert, and then normalize geodist maps
% invdist = distbroca.*broca;
% invdist = invdist*-1;
% normdist = (invdist - min(invdist)) / (max(invdist) - min(invdist));
% 
% new45 = results45.*normdist';
% new44 = results44.*normdist';

%% apply spatial constraint by weighting partial correlation results by manual probability maps within vlpfc

% prob44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_mean_101.1D']);
% prob45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_mean_101.1D']);
% 
% min44 = min(prob44(prob44>0));
% prob44 = prob44(find(broca));
% prob44(prob44==0)=min44; % set zeros within broca to minimum nonzero value
% new44 = zeros(length(M),1);
% new44(find(broca)) = prob44;
% normprob44 = (new44 - min(new44)) / (max(new44) - min(new44));
% normprob44 = normprob44*10;
% 
% min45 = min(prob45(prob45>0));
% prob45 = prob45(find(broca));
% prob45(prob45==0)=min45; % set zeros within broca to minimum nonzero value
% new45 = zeros(length(M),1);
% new45(find(broca)) = prob45;
% normprob45 = (new45 - min(new45)) / (max(new45) - min(new45));
% normprob45 = normprob45*10;
% 
% new45 = results45.*normprob45;
% new44 = results44.*normprob44;

% OR use prob map of both areas together to avoid gaps in transition zone
manprob = importdata(['/scr/murg2/MachineLearning/partialcorr/GroupProbBroca.1D']);
minprob = min(manprob(manprob>0));
manprob = manprob(find(broca));
manprob(manprob==0)=minprob; % set zeros within broca to minimum nonzero value
manprob = (log10(manprob) - min(log10(manprob))) ./ max(log10(manprob) - min(log10(manprob))); % change slope of prob map for less dramatic change
newprob = zeros(length(M),1);
newprob(find(broca)) = manprob;

new45 = results45.*newprob;
new44 = results44.*newprob;

% new45 = results45;
% new44 = results44;

% normnewprob = (newprob - min(newprob)) / (max(newprob) - min(newprob));
% normnewprob = normnewprob*10;
% new45 = results45.*normnewprob;
% new44 = results44.*normnewprob;

%% combine all partial correlation maps and apply winner-take-all

maps = [new45, new44, partcorrIC];
% maps = [results45, results44, partcorrIC];
% % normalize across nodes for each column
% maps_norm = zeros(size(maps));
% for i = 1:size(maps,2)
%     col = maps(:,i);
%     maps_norm(:,i) = (col - min(col)) / (max(col) - min(col));
% end
maps_norm = normc(maps); % normc preserves relative magnitude of columns by normalizing to a length of 1 instead of 0-1
% maps_norm=maps;

% apply winner-take-all to create binary partition from the partial correlation maps
[val,part] = max(maps_norm,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% save winner-take-all partition
filename = ['/scr/murg2/MachineLearning/partialcorr/' subject '_partcorrAll.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', part);
fclose(fid);

% keep only labels for BA44 and 45
part(part>2) = 0;

%% apply spatial constraint using cluster size (keep only largest clusters)

% get neighborhood information
surf = SurfStatReadSurf1(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/lh.very_inflated']);
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