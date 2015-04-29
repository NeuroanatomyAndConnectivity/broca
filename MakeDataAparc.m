% prepares HCP data for Broca classifier:
% extracts connectivity data from freesurfer opercularis + triangularis aparc labels and splits into three clusters: 0=neither, 1=BA45, 2=BA44
% includes additional rows of op and tri anatomical labels
% saves features and labels as separate variables

function[features, labels, op, tri] = MakeDataAparc(subject)

subject = num2str(subject);

%import freesurfer labels
AnatLabels = gifti(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/' subject '.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
vlpfc = op + tri;

% %import Broca labels
BA44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_' subject '.1D']);
BA45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_' subject '.1D']);

% %create label for "neither area"
neither = vlpfc + BA44 +BA45;
neither(neither>1)=0;

%read in the whole connectivity matrix as variable 'M'
fid = fopen(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% mask the full connectivity matrix by the manual labels
M_neither = M(find(neither),:);
M_BA45 = M(find(BA45),:);
M_BA44 = M(find(BA44),:);
% concatenate and save as features
features = horzcat(M_neither', M_BA45', M_BA44');

% create vectors of freesurfer labels:
op = op([find(neither); find(BA45); find(BA44)]);
tri = tri([find(neither); find(BA45); find(BA44)]);
op = op';
tri = tri';

% create vectors of the manual labels, where 0=neither, 1=BA45, 2=BA44
[R,C] = size(M_neither);
Vneither = zeros(1,R);
[R,C] = size(M_BA45);
V45 = ones(1,R);
[R,C] = size(M_BA44);
V44 = ones(1,R);
V44 = V44+1;
% concatenate and save as labels
labels = horzcat(Vneither, V45, V44);

