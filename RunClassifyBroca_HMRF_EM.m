sublist = importdata(['/scr/murg2/HCP500_glyphsets/subject_list.txt']);

radii = [];

for i=1:length(sublist)
radius = ClassifyBroca_HMRF_EM(sublist(i));

radii = horzcat(radii, radius);

fprintf('Subject # %u\n',i)

end

subradii = horzcat(sublist,radii')