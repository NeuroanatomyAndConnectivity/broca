function geoDist = surfGeoDist (surf, mask)
% Compute surface-based geodesic distance map from a given (gyral) mask 
% iterative backward algorithm that propagates along mesh topology 
% and that computes distance at each vertex v to the mask 
% if it has neighbor vertices for which the distance to the origin 
% has been computed in the previous iteration 
% 
% the distance is computed as the Euclidean distance along the edges of the mesh
% 
% > arguments: 
% inputmesh:    surf file as loaded by surfstat, containing surf.tri 
%               and surf.coord, and optionally surf.nbr
% mask:         can be any mask, 1 x num_vertices
%               vertices in the mask have 0 level depth
% > outputs: 
% geoDepth:     geodesic depth field from mask, 1 x num_vertices

% date:     Oct 11 2008
% authors:  boris@bic.mni.mcgill.ca
%           khs001@bic.mni.mcgill.ca
%

if (~isfield(surf,'nbr'))
    % requires surfGetNeighbors from MatlabTools/boris
    disp('Compute Neighbors using surfGetNeighbors.m')
    surf = surfGetNeighbors(surf);
    surf.nbr = surf.nbr'
else
    disp('Neighborhood info is already stored in surf.nbr');
end

origin = find(mask==1);

geoDepth         = zeros(1,length(surf.coord));

% indexPreCal(i)= 0 if the computation is not done yet on i-th vertex in
% previous iteration, indexPreCal(i)= 0 if computation is done. 
indexPreCal    = zeros(1,length(surf.coord));


% start the firefront from these vertices 
current_set    = origin;
% index the initial set
indexPreCal(current_set)=1;


k = 1;
k_1=0;
tmp=1;
% iteration #
l=2;
current_cal=surf.nbr(:,current_set);

% stop iteration when no computation is done in the previous step, 
% i.e., when k=k_1
disp('start iteration')
while k > k_1
    disp(['Iteration ' num2str(k)])
    k_1=k;
    
    for i=1:length(current_cal)
        t=current_cal(i);
        
        if (indexPreCal(t)==0)
            tmp=1000000;
            % if there is a neighbor with pre-computation
            
            for j=1:6
            
               % for all existing neigbors (sometimes there is 0 vertex nbr) 
               if surf.nbr(j,t) ~= 0
                   % compute depth if a current vertex has neighbor 
                   % vertices on which depth is computed in the 
                   % previous iteration 
                   
                   if ( indexPreCal(surf.nbr(j,t))>0 && ...
                        indexPreCal(surf.nbr(j,t))<l )
                       nbrCoord=surf.coord(:,surf.nbr(j,t));
                       curCoord=surf.coord(:,t);
                       dist_cur=sqrt(sum((nbrCoord - curCoord).^2,1));
                       geoDepth(t)=min(dist_cur+geoDepth(surf.nbr(j,t)), tmp);
                       tmp=geoDepth(t);
                       indexPreCal(t)=l;
                       % inculde neiborhood of the current vertex in
                       % next iteration 
                       current_set=[current_set surf.nbr(:,t)'];
                       
                   end
               end   
            end
           k=k+1;
        end   
    end
    
    % avoid double computation
    current_cal=unique(current_set);
    % remove 0 neighbor
    a=current_cal>0;
    current_cal=current_cal(a);
    current_set=[];
    l=l+1;
    %disp(l)
end

geoDist = geoDepth;