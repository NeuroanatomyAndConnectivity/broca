function[] = BrocaPartCorr(subjectID, datadir, outdir, varargin)

% Script to parcellate Broca's region into constituent areas 44 and 45 based on functional connetivity
% As described in:
%   Jakobsen, E., Liem, F., Klados, M., Bayrak, S., Petrides,
%   M., Margulies, D.S. Automated individual-level parcellation of
%   Broca's region based on functional connectivity. NeuroImage 2016
%   (online published)

% REQUIRED INPUTS
% 'subjectID' - subject identifier (e.g. 100307 for HCP data)
% 'datadir' - path to the directory where the input data are stored
    % should include:
        % individual subject directories named as subjectIDs containing:
            % correlation matrix file - rfMRI_REST_left_corr_avg.gii.data
            % anatomical freesurfer label file - <subjectID>.L.aparc.32k_fs_LR.label.gii
            % individual-level ICA results - ica_output_<subjectID>_20.mat
        % group connectivity maps - meanconn_nothresh_44.1D, meanconn_nothresh_45.1D
        % group probability maps - 44_mean_101.1D, 45_mean_101.1D
        % group-level ICA results - ica_output_20.mat
% 'outdir' - path to the desired output directory

% OPTIONAL ARGUMENTS
% 'rm_method' - method used to remove independent components. Can be:
    % 1 - (default) remove IC components by correlation threshold
        % 'corrthresh' - threshold used to remove components (default is r=0.4)
    % 2 - remove IC components by rank ordering (removes two most highly correlated components for each area)
% 'max_method' - method used to extract the individual connectivity templates from partial correlation results. Can be:
    % 1 - (default) use the single maximum partial correlation nodes
    % 2 - take the average of the top 5 percent maximum partial correlation nodes

% EXAMPLE
% BrocaPartCorr(100307, '/scr/murg2/HCP_Q3_glyphsets_left-only/', '/scr/murg2/MachineLearning/partialcorr/test/', 'corrthresh', 0.1, 'max_method', 2)

subjectID = num2str(subjectID);

% parse optional variables
i=1;
while i <= length(varargin)
  switch lower(varargin{i})
    case 'rm_method'
      rm_method = varargin{i+1}; i=i+2;
    case 'corrthresh'
      corrthresh = varargin{i+1}; i=i+2;
    case 'max_method'
       max_method = varargin{i+1}; i=i+2;
    otherwise
      error('Unknown option: %s\n',varargin{i}); i=i+1;
  end
end

% define default values
if ~exist('rm_method', 'var'), rm_method = 1; end
if ~exist('corrthresh', 'var'), corrthresh = 0.4; end
if ~exist('max_method', 'var'), max_method = 1; end

%% Import and mask data
disp('importing data...')

% get correlation matrix for individual subject
fid = fopen([datadir subjectID '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% get group connectivity maps
groupconn44 = importdata([datadir 'meanconn_nothresh_44.1D']);
groupconn44(isnan(groupconn44))=0;
groupconn45 = importdata([datadir 'meanconn_nothresh_45.1D']);
groupconn45(isnan(groupconn45))=0;

disp('masking data...')
% get pars opercularis and pars triangularis freesurfer labels for new dataset
AnatLabels = gifti([datadir subjectID '/' subjectID '.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
op = op';
tri = tri';
op = double(op);
tri = double(tri);

% import manual probability maps and binarize them
prob44 = importdata([datadir '44_mean_101.1D']);
prob45 = importdata([datadir '45_mean_101.1D']);
prob44_bin = prob44;
prob44_bin(prob44_bin>0)=1;
prob45_bin = prob45;
prob45_bin(prob45_bin>0)=1;

% add freesurfer labels to manual prob maps to create broca mask
broca = prob44_bin' + prob45_bin' + op + tri;
broca(broca>0)=1;
broca = logical(broca);

% mask the correlation matrix to Broca
M_small = M(:,broca);
M_small(isnan(M_small))=0;

%% Group-level partial correlation
disp('running partial correlations...')

% import group-level ICA maps
ica = load([datadir 'ica_output_20.mat']);
ica = ica.ic;

% run spatial correlation between group connectivity maps and ICA components
corr44 = zeros(0,size(ica,2));
corr45 = zeros(0,size(ica,2));
for i=1:size(ica,2)
    corr = corr2(groupconn44, ica(:,i));
    corr44(1,i) = corr;
    corr = corr2(groupconn45, ica(:,i));
    corr45(1,i) = corr;
end


if rm_method == 1 % remove components by correlation threshold
    rm = [find(corr45>corrthresh), find(corr44>corrthresh)];
    % rm = [find(abs(corr45)>corrthresh), find(abs(corr44)>corrthresh)];  % alternatively, take absolute value of r
elseif rm_method == 2 % OR remove components by correlation ranking
    [~, sort45] = sort(corr45, 'descend');
    [~, sort44] = sort(corr44, 'descend');
    rm = [sort45(:,1), sort45(:,2), sort44(:,1), sort44(:,2)];
end

ica_rm = ica;
ica_rm(:,rm) = [];

% also remove opposite group connectivity map
ica_rm45 = horzcat(ica_rm, groupconn44);
ica_rm44 = horzcat(ica_rm, groupconn45);

% run partial correlations with group connectivity maps
partcorr45 = partialcorr(M_small, groupconn45, ica_rm45);
partcorr44 = partialcorr(M_small, groupconn44, ica_rm44);

%% Extract individual-level connectivity templates
disp('extracting individual template maps...')

if max_method == 1 % extract the nodes with the maximum partial correlations to use as individual template maps
    [~,max45] = max(partcorr45);
    [~,max44] = max(partcorr44);
    indivconn45 = M_small(:,max45);
    indivconn44 = M_small(:,max44);
elseif max_method == 2 % OR average across the top 5 percent nodes
    topthresh45 = prctile(partcorr45(:),95);
    topnodes45 = partcorr45>topthresh45;
    topthresh44 = prctile(partcorr44(:),95);
    topnodes44 = partcorr44>topthresh44;
    indivconn45 = M_small(:,topnodes45);
    indivconn45 = mean(indivconn45,2);
    indivconn44 = M_small(:,topnodes44);
    indivconn44 = mean(indivconn44,2);
end

%% Individual-level partial correlation
disp('re-running partial correlations with individual template maps...')

% import individual-level ICA maps
ica_indiv = load([datadir subjectID '/ica_output_' subjectID '_20.mat']);
ica_indiv = ica_indiv.ic;

% run spatial correlation between individual connectivity maps and individual ICA components
corr44_indiv = zeros(0,size(ica_indiv,2));
corr45_indiv = zeros(0,size(ica_indiv,2));
for i=1:size(ica_indiv,2)
    corr = corr2(indivconn44, ica_indiv(:,i));
    corr44_indiv(1,i) = corr;
    corr = corr2(indivconn45, ica_indiv(:,i));
    corr45_indiv(1,i) = corr;
end

if rm_method == 1 % remove components by correlation threshold
    rm = [find(corr45_indiv>corrthresh), find(corr44_indiv>corrthresh)];
    % rm = [find(abs(corr45)>corrthresh), find(abs(corr44)>corrthresh)];  % alternatively, take absolute value of r
elseif rm_method == 2 % OR remove components by correlation ranking
    [~, sort45] = sort(corr45_indiv, 'descend');
    [~, sort44] = sort(corr44_indiv, 'descend');
    rm = [sort45(:,1), sort45(:,2), sort44(:,1), sort44(:,2)];
end

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
partcorrIC = zeros(size(M,1), size(ica_rm_indiv,2));
for i=1:size(ica_rm_indiv,2)
    rm = ica_rm_indiv;
    rm(:,i) = [];
    rm = horzcat(indivconn45, indivconn44, rm);
    partcorr = partialcorr(M_small, ica_rm_indiv(:,i), rm);
    results = zeros(length(M),1);
    results(broca) = partcorr;
    partcorrIC(:,i) = results;
end
partcorrIC(isnan(partcorrIC))=0;

% make results wholebrain
results45 = zeros(length(M),1);
results45(broca) = partcorr45;
results44 = zeros(length(M),1);
results44(broca) = partcorr44;

%% Spatial weighting
disp('applying spatial weighting...')

% weight partial correlation results by manual probability maps within vlpfc
min44 = min(prob44(prob44>0));
prob44 = prob44(broca);
prob44(prob44==0)=min44; % set zeros within broca to minimum nonzero value
prob44 = (log10(prob44) - min(log10(prob44))) ./ (max(log10(prob44)) - min(log10(prob44)));
norm44 = zeros(length(M),1);
norm44(broca) = prob44;

min45 = min(prob45(prob45>0));
prob45 = prob45(broca);
prob45(prob45==0)=min45; % set zeros within broca to minimum nonzero value
prob45 = (log10(prob45) - min(log10(prob45))) ./ (max(log10(prob45)) - min(log10(prob45)));
norm45 = zeros(length(M),1);
norm45(broca) = prob45;

new45 = results45.*norm45;
new44 = results44.*norm44;

%% Winner-take-all partition
disp('applying winner-take-all partition...')

%combine all partial correlation maps
maps = [new45, new44, partcorrIC];

% apply winner-take-all to create binary partition from the partial correlation maps
[val,part] = max(maps,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% save winner-take-all partition
filename = sprintf('%s%s_BrocaPart_WTA.',outdir, subjectID);
fid = fopen([filename '1D'],'w');
fprintf(fid, '%u\n', part);
fclose(fid);

%% Extract largest clusters labeled as 44 and 45
disp('extracting largest clusters...')

% keep only labels for BA44 and 45
part(part>2) = 0;

% get neighborhood information
surf = SurfStatReadSurf1([datadir subjectID '/lh.very_inflated']);
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
disp('saving results...')

% combine results from BA 44 and 45
results = results45 + results44;
filename = sprintf('%s%s_BrocaPart.',outdir, subjectID);
fid = fopen([filename '1D'],'w');
fprintf(fid, '%u\n', results);
fclose(fid);

figure('visible','off'); SurfStatViewDataLat(results,surf);
saveas(gcf, [filename 'png']); close all;

disp('done.')

