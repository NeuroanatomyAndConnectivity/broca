%%  The MAP algorithm
%---input---------------------------------------------------------
%   X: initial 2D labels
%   Y: image
%   Z: 2D constraints
%   mu: vector of means
%   sigma: vector of standard deviations
%   k: number of labels
%   MAP_iter: maximum number of iterations of the MAP algorithm
%   show_plot: 1 for showing a plot of energy in each iteration
%       and 0 for not showing
%---output--------------------------------------------------------
%   X: final 2D labels
%   sum_U: final energy

%   Copyright by Quan Wang, 2012/04/25
%   Please cite: Quan Wang. HMRF-EM-image: Implementation of the 
%   Hidden Markov Random Field Model and its Expectation-Maximization 
%   Algorithm. arXiv:1207.3510 [cs.CV], 2012.

function [X sum_U]=MRF_MAP_surf(X,Y,Z,mu,sigma,k,MAP_iter,show_plot,edg, distances, radius)

x=X(:);
y=Y(:);
m=length(Y);
U=zeros(m,k);
sum_U_MAP=zeros(1,MAP_iter);
for ot=1:MAP_iter % iterations
    fprintf('  Inner iteration: %d\n',ot);
    U1=U;
    U2=U;
    
    for l=1:k % all labels
        yi=y-mu(l);
        temp1=yi.*yi/sigma(l)^2/2;
        temp1=temp1+log(sigma(l));
        U1(:,l)=U1(:,l)+temp1;	        

        for ind=1:m % all nodes    
		    u2=0;
		    % edg could be replaced with the node indices within r = radius (mm).
		    % e = unique(nonzeros(edg(:,nonzeros(edg(:,ind)))));
		    %if l == 3
            %    e = find(distances(:,ind) < radius/2);
            %else
                e = find(distances(:,ind) < radius);
            %end
		    for j = 1:length(e)		    	
		        if Z(e(j))==0
                    % probability values could be used here
		            u2=u2+(l ~= X(e(j)))/2;
		        end	
		    end    
		    U2(ind,l)=u2;
		end
    end
    U=U1+U2;
    [temp x]=min(U,[],2);
    sum_U_MAP(ot)=sum(temp(:));
    
    X=x(:);
    if ot>=3 && std(sum_U_MAP(ot-2:ot))/sum_U_MAP(ot)<0.0001
        break;
    end
end

sum_U=0;
for ind=1:m % all pixels
    sum_U=sum_U+U(ind,x(ind));
end
if show_plot==1
    figure;
    plot(1:ot,sum_U_MAP(1:ot),'r');
    title('sum U MAP');
    xlabel('MAP iteration');
    ylabel('sum U MAP');
    drawnow;
end
