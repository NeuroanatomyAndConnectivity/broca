function[results] = BrocaPartCorrICA10_NKI(subject)

subject = char(subject);
subdir = ('/scr/murg2/MachineLearning/partialcorr/NKI/subjects/');
outdir = ('/scr/murg2/MachineLearning/partialcorr/NKI/results/');

%% Import data and run partial correlation for each area

% get correlation matrix for individual subject
fid = fopen([subdir subject '_lh_preprocessed_fsaverage5_fwhm6_corr.bin'], 'r');
M = fread(fid,[10242 10242], 'float32');

% get group connectivity maps
groupconn44 = gifti(['/scr/murg2/MachineLearning/partialcorr/NKI/BrocaGroupMaps_fsa5/groupconn44_fsa5.func.gii']);
groupconn44 = groupconn44.cdata;
groupconn44(isnan(groupconn44))=0;
groupconn45 = gifti(['/scr/murg2/MachineLearning/partialcorr/NKI/BrocaGroupMaps_fsa5/groupconn45_fsa5.func.gii']);
groupconn45 = groupconn45.cdata;
groupconn45(isnan(groupconn45))=0;

% get pars opercularis and pars triangularis freesurfer labels
[~, label, ~] = read_annotation('/scr/murg2/MachineLearning/partialcorr/NKI/lh.aparc.annot');
op = label == 9221340;
tri = label == 1326300;
op = op';
tri = tri';
op = double(op);
tri = double(tri);

% import manual probability maps and binarize them
prob44 = gifti(['/scr/murg2/MachineLearning/partialcorr/NKI/BrocaGroupMaps_fsa5/groupprob44_fsa5.func.gii']);
prob44 = prob44.cdata;
prob45 = gifti(['/scr/murg2/MachineLearning/partialcorr/NKI/BrocaGroupMaps_fsa5/groupprob45_fsa5.func.gii']);
prob45 = prob45.cdata;
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

% import ICA maps
ica = load('/scr/murg2/MachineLearning/partialcorr/NKI/ICA10_fsa5/ICA10_fsa5.mat');
ica = ica.ICA10_fsa5;

ica = load('/scr/murg2/MachineLearning/partialcorr/NKI/ICA10_NKI.mat');
ica = ica.A_ts;

% run spatial correlation between group connectivity maps and ICA components
corr44 = [];
corr45 = [];
for i=1:size(ica,2)
    corr = corr2(groupconn44, ica(:,i));
    corr44 = horzcat(corr44, corr);
    corr = corr2(groupconn45, ica(:,i));
    corr45 = horzcat(corr45, corr);
end

% remove the two components most highly correlated with each group connectivity map
[~, sort45] = sort(corr45, 'descend');
[~, sort44] = sort(corr44, 'descend');
rm = [sort45(:,1), sort45(:,2), sort44(:,1), sort44(:,2)];

% % OR remove components with correlations over a certain value
% rm = [find(corr45>0.2), find(corr44>0.2)];

ica_rm = ica;
ica_rm(:,rm) = [];

% also remove opposite group connectivity map
ica_rm45 = horzcat(ica_rm, groupconn44);
ica_rm44 = horzcat(ica_rm, groupconn45);

% run partial correlations with group connectivity maps
partcorr45 = partialcorr(M_small, groupconn45, ica_rm45);
partcorr44 = partialcorr(M_small, groupconn44, ica_rm44);

% find index of max nodes from partcorr maps i.e. nodes with most similar connectivity to group-level BA44 and 45 maps
[~,max45] = max(partcorr45);
[~,max44] = max(partcorr44);

% extract connectivity maps from max nodes, then use as new template maps and re-run partial correlation
indivconn45 = M_small(:,max45);
indivconn44 = M_small(:,max44);
ica_rm45 = horzcat(ica_rm, indivconn44);
ica_rm44 = horzcat(ica_rm, indivconn45);
partcorr45 = partialcorr(M_small, indivconn45, ica_rm45);
partcorr44 = partialcorr(M_small, indivconn44, ica_rm44);
partcorr45(isnan(partcorr45))=0;
partcorr44(isnan(partcorr44))=0;

%run partial correlations for remaining ICA components
partcorrIC = [];
for i=1:size(ica_rm,2)
    rm = ica_rm;
    rm(:,i) = [];
    rm = horzcat(groupconn45, groupconn44, rm);
    partcorr = partialcorr(M_small, ica_rm(:,i), rm);
    results = zeros(length(M),1);
    results(find(broca)) = partcorr;
    partcorrIC = horzcat(partcorrIC, results);
end

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

%% apply spatial constraint by weighting partial correlation results by manual probability maps within vlpfc

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

%% combine all partial correlation maps and apply winner-take-all

maps = [new45, new44, partcorrIC];
maps_norm = normc(maps); % normc preserves relative magnitude of columns by normalizing to a length of 1 instead of 0-1

% apply winner-take-all to create binary partition from the partial correlation maps
[val,part] = max(maps_norm,[],2); %return the maximum from each row (returns 1 if values are the same)
part(val==0) = 0;

% save winner-take-all partition
filename = [outdir subject '_partcorrAll_rmtop2.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', part);
fclose(fid);

% keep only labels for BA44 and 45
part(part>2) = 0;

%% extract largest clusters

% get neighborhood information
surf = SurfStatReadSurf1('/scr/murg2/MachineLearning/partialcorr/NKI/lh.inflated');
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

% combine results from BA 44 and 45
results = results45 + results44;
filename = [outdir subject '_partcorrMaxclust_rmtop2.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);

figure('visible','off'); SurfStatView(results,surf);
saveas(gcf, [outdir subject '_partcorrMaxclust_rmtop2.png']); close all;