data = load('/SCR2/HCP500_glyphsets/results/results_991267.1D');
surf_gii = gifti(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii']);
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;

aparc = gifti(['/scr/dattel2/100307/MNINonLinear/fsaverage_LR32k/100307.L.aparc.a2009s.32k_fs_LR.label.gii']);
aparc = aparc.cdata;
roi = zeros(32492,1);
roi(find(aparc == 14 | aparc == 69 | aparc == 12)) = 1;
runIt = 1;
while runIt == 1
    label = data;
    clear clus
    edg = SurfStatEdg(surf);
    lab = nonzeros(unique(label));
    clus.label = zeros([1 length(surf.coord)]);
    countC = 0;
    clus.label = zeros([1 32492]);
    clus.network = zeros([1 32492]);
    for i = 1:length(lab)    
        a = zeros([length(surf.coord) 1]);
        a(find(label == lab(i))) = 1; 
        slm = struct();
        slm.tri = surf.tri';
        slm.t = a';
        [cluster,clusid] = SurfStatPeakClus(slm,ones([length(surf.coord) 1]),0.05, ones(1,length(surf.coord)), edg);  
        clear slm
        for j = 1:length(clusid.clusid)
            nodes = cluster.vertid(find(cluster.clusid == j));  
                countC = countC + 1;
                clus.label(nodes) = countC;   
                clus.network(nodes) = lab(i);       
                clus.nverts(nodes) = clusid.nverts(j);
        end
        clear cluster clusid
    end

    d1 = clus.nverts(clus.network == 1);
    max1 = max(d1(find(roi(clus.network == 1))));
    d2 = clus.nverts(clus.network == 2);
    max2 = max(d2(find(roi(clus.network == 2))));

    modClus = clus.network;
    data1 = zeros(32492,1);
    data1(find(clus.network == 1)) = clus.nverts(find(clus.network == 1));
    modClus(find(data1 ~= max1 & data1 ~= 0)) = 2;

    data2 = zeros(32492,1);
    data2(find(clus.network == 2)) = clus.nverts(find(clus.network == 2));
    modClus(find(data2 ~= max2 & data2 ~= 0)) = 1;

    if (length(find(data1 ~= max1 & data1 ~= 0)) + length(find(data2 ~= max2 & data2 ~= 0))) == 0
        runIt = 0;
        disp('Done!');
    else
        data = modClus;
        disp(['running again. ' num2str((length(find(data1 ~= max1 & data1 ~= 0)) + length(find(data2 ~= max2 & data2 ~= 0)))) ' nodes changed.']);
    end
end

modClusROI = modClus;
modClusROI(find(roi == 0)) = 0;


label = modClusROI;
clear clus
edg = SurfStatEdg(surf);
lab = nonzeros(unique(label));
clus.label = zeros([1 length(surf.coord)]);
countC = 0;
clus.label = zeros([1 32492]);
clus.network = zeros([1 32492]);
for i = 1:length(lab)    
    a = zeros([length(surf.coord) 1]);
    a(find(label == lab(i))) = 1; 
    slm = struct();
    slm.tri = surf.tri';
    slm.t = a';
    [cluster,clusid] = SurfStatPeakClus(slm,ones([length(surf.coord) 1]),0.05, ones(1,length(surf.coord)), edg);  
    clear slm
    for j = 1:length(clusid.clusid)
        nodes = cluster.vertid(find(cluster.clusid == j));  
            countC = countC + 1;
            clus.label(nodes) = countC;   
            clus.network(nodes) = lab(i);       
            clus.nverts(nodes) = clusid.nverts(j);
    end
    clear cluster clusid
end

d1 = clus.nverts(clus.network == 1);
max1 = max(d1(find(roi(clus.network == 1))));
d2 = clus.nverts(clus.network == 2);
max2 = max(d2(find(roi(clus.network == 2))));

modClusROI2 = clus.network;
data1 = zeros(32492,1);
data1(find(clus.network == 1)) = clus.nverts(find(clus.network == 1));
modClusROI2(find(data1 ~= max1 & data1 ~= 0)) = 0;

data2 = zeros(32492,1);
data2(find(clus.network == 2)) = clus.nverts(find(clus.network == 2));
modClusROI2(find(data2 ~= max2 & data2 ~= 0)) = 0;