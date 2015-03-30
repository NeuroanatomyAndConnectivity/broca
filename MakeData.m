function[data] = MakeData(subject)

subject = num2str(subject);

BA44 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/44_' subject '.1D']);
BA45 = importdata(['/scr/murg2/HCP_Q3_glyphsets_left-only/post-Montreal_labels/45_' subject '.1D']);

%read in the whole connectivity matrix as variable 'M'
fid = fopen(['/scr/murg2/HCP_Q3_glyphsets_left-only/' subject '/rfMRI_REST_left_corr_avg.gii.data'], 'r');
M = fread(fid,[32492 32492], 'float32');

% mask the full connectivity matrix by the manual labels and concatenate with additional row indicating labels
M_BA45 = M(find(BA45),:);
[R,C] = size(M_BA45);
V45 = ones(1,R);
data_BA45 = vertcat(M_BA45', V45);

M_BA44 = M(find(BA44),:);
[R,C] = size(M_BA44);
V44 = ones(1,R);
V44 = V44+1;
data_BA44 = vertcat(M_BA44', V44);

data = horzcat(data_BA45, data_BA44);