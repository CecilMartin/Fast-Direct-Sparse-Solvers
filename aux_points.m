% An auxiliary file to help illustrate the ordering of points
n = 10;
idx = zeros(n+1,n+1);  % index mapping to each point, including "ghost" points
idx(2:n,2:n) = reshape(1:(n-1)^2,n-1,n-1);
[x1, x2] = ndgrid((1:(n-1))/n); x = [x2(:) 1-x1(:)];

% HERE idx shows the indexing of the uniform grid 
% while x gives out the coordinate of the points under the given indexing