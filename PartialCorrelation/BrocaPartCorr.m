function[] = BrocaPartCorr(subject, corrthresh, fileext)

subject = num2str(subject);
subdir = ('/scr/murg2/HCP_Q3_glyphsets_left-only/');
% subdir = ('/scr/murg1/HCP500_glyphsets/');
outdir = ('/scr/murg2/MachineLearning/partialcorr/20comps_results/HCP_ICA_Indiv/');
% outdir = ('/scr/murg2/MachineLearning/partialcorr/20comps_results/HCP_ICA_Indiv/HCP150_new/');

%% Import and mask data

% get correlation matrix for individual subject
fid = fopen([subdir subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% get group connectivity maps
groupconn44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_nothresh_44.1D']);
groupconn44(isnan(groupconn44))=0;
groupconn45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_maps/meanconn_nothresh_45.1D']);
groupconn45(isnan(groupconn45))=0;

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
prob44_bin = prob44;
prob44_bin(prob44_bin>0)=1;
prob45_bin = prob45;
prob45_bin(prob45_bin>0)=1;

% add freesurfer labels to manual prob maps to create broca mask
broca = prob44_bin' + prob45_bin' + op + tri;
broca(broca>0)=1;

% mask the correlation matrix to Broca
M_small = M(:,find(broca));
M_small(isnan(M_small))=0;

%% Group-level partial correlation

% import group-level ICA maps
ica = load('/scr/murg2/MachineLearning/partialcorr/ICA/ICA_HCP/ica_output_20.mat');
ica = ica.ic;
% surf = SurfStatReadSurf1(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/lh.inflated']);

% run spatial correlation between group connectivity maps and ICA components
corr44 = [];
corr45 = [];
for i=1:size(ica,2)
    corr = corr2(groupconn44, ica(:,i));
    corr44 = horzcat(corr44, corr);
    corr = corr2(groupconn45, ica(:,i));
    corr45 = horzcat(corr45, corr);
end

% remove components by correlation ranking
% [~, sort45] = sort(corr45, 'descend');
% [~, sort44] = sort(corr44, 'descend');
% rm = [sort45(:,1), sort45(:,2), sort44(:,1), sort44(:,2)];

% % OR remove components by threshold
rm = [find(corr45>corrthresh), find(corr44>corrthresh)];
% taking absolute value of r
% rm = [find(abs(corr45)>corrthresh), find(abs(corr44)>corrthresh)];

ica_rm = ica;
ica_rm(:,rm) = [];

% also remove opposite group connectivity map
ica_rm45 = horzcat(ica_rm, groupconn44);
ica_rm44 = horzcat(ica_rm, groupconn45);

% run partial correlations with group connectivity maps
partcorr45 = partialcorr(M_small, groupconn45, ica_rm45);
partcorr44 = partialcorr(M_small, groupconn44, ica_rm44);

%% Extract individual-level connectivity templates

% find index of max nodes from partcorr maps i.e. nodes with most similar connectivity to group-level BA44 and 45 maps
[~,max45] = max(partcorr45);
[~,max44] = max(partcorr44);

% extract connectivity maps from max nodes to use as new template maps
indivconn45 = M_small(:,max45);
indivconn44 = M_small(:,max44);

% OR use the average of the top 5 percent as individual template
% %find indices of the top 5 percent partial correlations
% topthresh45 = prctile(partcorr45(:),95);
% topnodes45 = find(partcorr45>topthresh45);
% topthresh44 = prctile(partcorr44(:),95);
% topnodes44 = find(partcorr44>topthresh44);
% 
% % extract connectivity maps from max nodes to use as new template maps
% indivconn45 = M_small(:,topnodes45);
% indivconn45 = mean(indivconn45,2);
% indivconn44 = M_small(:,topnodes44);
% indivconn44 = mean(indivconn44,2);


%% Individual-level partial correlation

% import individual-level ICA maps
ica_indiv = load(['/scr/murg2/MachineLearning/partialcorr/ICA/ICA_HCP/ica_output_' subject '_20.mat']);
ica_indiv = ica_indiv.ic;

% run spatial correlation between individual connectivity maps and individual ICA components
corr44_indiv = [];
corr45_indiv = [];
for i=1:size(ica_indiv,2)
    corr = corr2(indivconn44, ica_indiv(:,i));
    corr44_indiv = horzcat(corr44_indiv, corr);
    corr = corr2(indivconn45, ica_indiv(:,i));
    corr45_indiv = horzcat(corr45_indiv, corr);
end

% % remove components by correlation ranking
% [~, sort45] = sort(corr45_indiv, 'descend');
% [~, sort44] = sort(corr44_indiv, 'descend');
% rm = [sort45(:,1), sort45(:,2), sort44(:,1), sort44(:,2)];

% % OR remove components by threshold
rm = [find(corr45_indiv>corrthresh), find(corr44_indiv>corrthresh)];
% taking absolute value of r
% rm = [find(abs(corr45_indiv)>corrthresh), find(abs(corr44_indiv)>corrthresh)];

ica_rm_indiv = ica_indiv;
ica_rm_indiv(:,rm) = [];

% also regress out opposite group connectivity map
ica_rm45_indiv = horzcat(ica_rm_indiv, indivconn44);
ica_rm44_indiv = horzcat(ica_rm_indiv, indivconn45);

% re-run partial correlation with individual-level template maps
partcorr45 = partialcorr(M_small, indivconn45, ica_rm45_indiv);
partcorr44 = partialcorr(M_small, indivconn44, ica_rm44_indiv);
partcorr45(isnan(partcorr45))=0;
partcorr44(isnan(partcorr44))=0;

%run partial correlations for remaining ICA components
partcorrIC = [];
for i=1:size(ica_rm_indiv,2)
    rm = ica_rm_indiv;
    rm(:,i) = [];
    rm = horzcat(indivconn45, indivconn44, rm);
    partcorr = partialcorr(M_small, ica_rm_indiv(:,i), rm);
    results = zeros(length(M),1);
    results(find(broca)) = partcorr;
    partcorrIC = horzcat(partcorrIC, results);
end
partcorrIC(isnan(partcorrIC))=0;

% make results wholebrain and save as 1D
results45 = zeros(length(M),1);
results45(find(broca)) = partcorr45;
% filename = [outdir subject '_partcorr45_indiv.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', results45);
% fclose(fid);

results44 = zeros(length(M),1);
results44(find(broca)) = partcorr44;
% filename = [outdir subject '_partcorr44_indiv.1D'];
% fid = fopen(filename,'w');
% fprintf(fid, '%u\n', results44);
% fclose(fid);

%% Spatial weighting

% weight partial correlation results by manual probability maps within vlpfc
min44 = min(prob44(prob44>0));
prob44 = prob44(find(broca));
prob44(prob44==0)=min44; % set zeros within broca to minimum nonzero value
prob44 = (log10(prob44) - min(log10(prob44))) ./ (max(log10(prob44)) - min(log10(prob44)));
norm44 = zeros(length(M),1);
norm44(find(broca)) = prob44;

min45 = min(prob45(prob45>0));
prob45 = prob45(find(broca));
prob45(prob45==0)=min45; % set zeros within broca to minimum nonzero value
prob45 = (log10(prob45) - min(log10(prob45))) ./ (max(log10(prob45)) - min(log10(prob45)));
norm45 = zeros(length(M),1);
norm45(find(broca)) = prob45;

new45 = results45.*norm45;
new44 = results44.*norm44;

%% Winner-take-all partition

%combine all partial correlation maps
maps = [new45, new44, partcorrIC];

% apply winner-take-all to create binary partition from the partial correlation maps
[val,part] = max(maps,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% save winner-take-all partition
filename = [outdir subject '_ICA_indiv_WTA_SW_' fileext '.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', part);
fclose(fid);

%% Extract largest clusters labeled as 44 and 45

% keep only labels for BA44 and 45
part(part>2) = 0;

% get neighborhood information
surf = SurfStatReadSurf1([subdir subject '/lh.very_inflated']);
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
[~,~,members45] = networkComponents(A45);
% find the largest component of BA45
largest45 = members45{1};
results45 = zeros(size(part));
results45(largest45) = 1;

% create edge list for BA44
label44 = part';
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
[~,~,members44] = networkComponents(A44);
% find the largest component of BA44
largest44 = members44{1};
results44 = zeros(size(part));
results44(largest44) = 2;

%% Save results

% combine results from BA 44 and 45
results = results45 + results44;
filename = [outdir subject '_ICA_indiv_SW_' fileext '.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);

figure('visible','off'); SurfStatView(results,surf);
saveas(gcf, [outdir subject '_ICA_indiv_SW_' fileext '.png']); close all;

