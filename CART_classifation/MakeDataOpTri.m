% prepares HCP data for Broca classifier:
% extracts connectivity data from freesurfer opercularis + triangularis aparc labels and splits into three clusters: 0=neither, 1=BA45, 2=BA44
% includes additional rows of op and tri anatomical labels
% saves features and labels as separate variables

function[op, tri] = MakeDataOpTri(subject)

subject = num2str(subject);

%import freesurfer labels
AnatLabels = gifti(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/' subject '.L.aparc.32k_fs_LR.label.gii']);
AnatLabelsData = AnatLabels.cdata;
op = AnatLabelsData == 18;
tri = AnatLabelsData == 20;
vlpfc = op + tri;
% get surface for new dataset
%may need to convert surface from .asc, use freesurfer mris_convert L.midthickness.asc lh.midthickness
surf = SurfStatReadSurf1(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/lh.midthickness']);
surf = surfGetNeighbors(surf);

% calculate geodesic distance from freesurfer labels
distTri = surfGeoDist_parcellation(surf, tri');
distOp = surfGeoDist_parcellation(surf, op');

% %import Broca labels
BA44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_' subject '.1D']);
BA45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_' subject '.1D']);

% %create label for "neither area"
neither = vlpfc + BA44 +BA45;
neither(neither>1)=0;

% create vectors of freesurfer labels:
op = distOp([find(neither); find(BA45); find(BA44)]);
tri = distTri([find(neither); find(BA45); find(BA44)]);
op = op';
tri = tri';
% % add these labels to features variable
% features = vertcat(op', tri', features);


