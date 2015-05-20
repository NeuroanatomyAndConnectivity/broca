function surf = surfGetNeighbors(surf)
% surf = surfGetNeighbors(surf)
% compute neigbors of each vertex in a given surf file
% > arguments:
%               surf - a surface struct containing surf.tri and surf.coord
% > outputs: 
%               surf -  a surface struct containing surf.tri, surf.coord, 
%                       and surf.nbr

% taken from: Prof. Moo Chungs website, with slight modifications
% http://www.stat.wisc.edu/~mchung/
% boris@bic.mni.mcgill.ca


% useful variables
surf.tri = surf.tri';
num_tri = size(surf.tri,2);
num_points = length(surf.coord);


% compute the maximum degree of node
degree=zeros(1,num_points);
for j=1:num_tri
    degree(surf.tri(:,j)) = degree(surf.tri(:,j))+1;
end
max_degree=max(degree);


% set up neigboring nodes field in surf
nbr=zeros(max_degree,num_points);
for i_tri=1:num_tri
    for j=1:3
        cur_point = surf.tri(j,i_tri);
        for k=1:3
            if (j ~= k)
                nbr_point= surf.tri(k,i_tri);
                if find(nbr(:,cur_point)==nbr_point)
                    ;
                else
                    n_nbr = min(find(nbr(:,cur_point) == 0));
                    nbr(n_nbr,cur_point) = nbr_point;
                end;
            end;
        end;
    end;
end;

surf.tri = surf.tri';
% matrix with neighbor information 
% surf.nbr(:,i): means 6 adjacent vertices of i-th vertex
surf.nbr = nbr;