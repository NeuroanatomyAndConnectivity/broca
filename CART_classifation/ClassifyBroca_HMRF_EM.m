% sublist = importdata(['/scr/murg2/HCP500_glyphsets/subject_list.txt']);
% for i=1:length(sublist)
% ClassifyBroca_HMRF_EM(sublist(i));
% end


function[radius Xa X mu sigma surf] = ClassifyBroca_HMRF_EM(subject)

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

% mask features
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

% make results same size as surface and save
results = zeros([length(broca) 1]);
results(find(broca)) = label;

filename = ['/scr/murg2/MachineLearning/HMRF/results/' subject '_orig.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', results);
fclose(fid);

surf = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
figure('visible','off'); SurfStatView(results,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/HMRF/figures/' subject '_orig_wholebrain.png']); close all;

% make probability maps for each area
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

BA45prob = zeros([length(broca) 1]);
BA45prob(find(broca)) = probmap45;

Neitherprob = zeros([length(broca) 1]);
Neitherprob(find(broca)) = probmapNeither;

%% Apply HMRF_EM framework to improve spatial continuity of labels

%   Modified from: Quan Wang. HMRF-EM-image: Implementation of the 
%   Hidden Markov Random Field Model and its Expectation-Maximization 
%   Algorithm. arXiv:1207.3510 [cs.CV], 2012.

%---input---------------------------------------------------------
%   X: initial binarized labels
%   Y: probability values (rows = indicies; cols = prob vals)
%   Z: hard binarized edges
%   mu: initial vector of means
%   sigma: initial vector of standard deviations
%   k: 3
%   EM_iter: maximum number of iterations of the EM algorithm
%   MAP_iter: maximum number of iterations of the MAP algorithm
%---output--------------------------------------------------------
%   X: final 2D labels
%   mu: final vector of means
%   sigma: final vector of standard deviations

%% Set parameters
% [broca] = find(results>0);
% brocasize = size(broca);
% radius = round((brocasize(1,1))*0.01); % search radius (neighborhood size) in number of edges currently set to 1 percent of size of broca
radius = 4; % search radius (neighborhood size) in number of edges
EM_iter = 10; % maximum number of iterations for EM algorithm
MAP_iter = 10; % maximum number of iterations for MAP algorithm
k = 3; % number of classes

%% Get input data
surf = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);

% Define X: initial binarized labels
%X = load(['/scr/murg2/MachineLearning/newdata_results/' subject '_class_prob_nbrs.1D']); 
X = results;

X(X == 0) = 3; % set neither to 3

% Define Z: hard binarized edges
Z = zeros(length(X),1);

% Define Y: probability values (rows = indices; columns = prob values)
% a = load(['/scr/murg2/MachineLearning/newdata_results/' subject '_class_probmap45.1D']);
% b = load(['/scr/murg2/MachineLearning/newdata_results/' subject '_class_probmap44.1D']);
% c = load(['/scr/murg2/MachineLearning/newdata_results/' subject '_class_probmapNeither.1D']);
a = BA45prob;
b = BA44prob;
c = Neitherprob;
c = c/2; % !!!reducing probability of neither!!!
Y = [a b c];
Y = reshape(Y,32492,1,3);
Y = rgb2gray(Y);
roi = find(Y);

% unsure of whether this should be done! Too much of area 6 is included with this step
% a = SurfStatSmooth(a',surf, 1)';
% b = SurfStatSmooth(b',surf, 1)';
% c = SurfStatSmooth(c',surf, 1)';
% [~,X] = max([a b c],[],2); 
% X(find(a+b == 0)) = 3;
% X(find(Y == 0)) = 0;

Y = Y(roi);
X = X(roi);
Z = Z(roi);

surf.coord = surf.coord(:,roi);
surf.tri = surf.tri(find(ismember(surf.tri(:,1),roi)+ismember(surf.tri(:,2),roi)+ismember(surf.tri(:,3),roi) == 3),:);

for i = 1:length(roi); 
    surf.tri(find(surf.tri(:,1) == roi(i)),1) = i; 
    surf.tri(find(surf.tri(:,2) == roi(i)),2) = i; 
    surf.tri(find(surf.tri(:,3) == roi(i)),3) = i; 
end

% get neighbors
edg = surfGetNeighbors(surf);

% visualize and save image of uncorrected labels
figure('visible','off'); SurfStatView(X,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/HMRF/figures/' subject '_orig.png']); close all;
distances = squareform(pdist(surf.coord'));

% calculate mu (initial vector of means) and sigma (initial vector of standard deviations)
y=Y(:);
mu=zeros(k,1);
sigma=zeros(k,1);
for i=1:k
    yy=y(X==i);
    mu(i)=mean(yy);
    sigma(i)=std(yy);
end

%% Begin iterations:

m=length(Y);
P_lyi=zeros(k,m);
sum_U=zeros(1,EM_iter);

for it=1:EM_iter
    fprintf('Iteration: %d\n',it);
% MAP: 
    [X, sum_U(it)]=MRF_MAP_surf(X,Y,Z,mu,sigma,k,MAP_iter,0, edg, distances, radius);	    
    x=X(:);
    %update mu and sigma    
    % get P_lyi
    for l=1:k % all labels
        temp1=1/sqrt(2*pi*sigma(l)^2)*exp(-(y-mu(l)).^2/2/sigma(l)^2);
        temp2=temp1*0;
        for ind=1:m % all nodes    
            u=0;
            % edg could be replaced with the node indices within r = radius (mm).
            % currently grabbing 3 edge radius:		    
            % e = unique(nonzeros(edg(:,nonzeros(edg(:,ind)))));
            %if l == 3
            %    e = find(distances(:,ind) < radius/2);
            %else
                e = find(distances(:,ind) < radius);
            %end
            for j = 1:length(e)		    	
                if Z(e(j))==0
                    % probability values could be used here
                    u=u+(l ~= X(e(j)))/2;
                end		    
            end    
            temp2(ind)=u;
        end

        P_lyi(l,:)=temp1.*exp(-temp2);
    end
    temp3=sum(P_lyi,1);
    P_lyi=bsxfun(@rdivide,P_lyi,temp3);

    % get mu and sigma
    for l=1:k % all labels
        mu(l)=P_lyi(l,:)*y;
        mu(l)=mu(l)/sum(P_lyi(l,:));
        sigma(l)=P_lyi(l,:) * ( (y-mu(l)).^2 );
        sigma(l)=sigma(l)/sum(P_lyi(l,:));
        sigma(l)=sqrt(sigma(l));
    end


    if it>=3 && std(sum_U(it-2:it))/sum_U(it)<0.0001
        break;
    end
end

%% Visualize and save results
figure('visible','off');    
plot(1:it,sum_U(1:it),'LineWidth',2);
hold on;
plot(1:it,sum_U(1:it),'.','MarkerSize',20);
title('sum of U in each EM iteration');
xlabel('EM iteration');
ylabel('sum of U');
saveas(gcf, ['/scr/murg2/MachineLearning/HMRF/figures/' subject '_opimization.png']); close all;

% viz patch:
figure('visible','off');  SurfStatView(X,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/HMRF/figures/' subject '_postcorrection.png']); close all;
% viz wholebrain:
surf = SurfStatReadSurf1(['/scr/murg2/HCP500_glyphsets/' subject '/lh.very_inflated']);
Xa = zeros(length(surf.coord),1);
Xa(roi) = X;
figure('visible','off'); SurfStatView(Xa,surf);
saveas(gcf, ['/scr/murg2/MachineLearning/HMRF/figures/' subject '_wholebrain.png']); close all;

%save corrected labels as 1D
filename = ['/scr/murg2/MachineLearning/HMRF/results/' subject '_postcorrection.1D'];
fid = fopen(filename,'w');
fprintf(fid, '%u\n', Xa);
fclose(fid);