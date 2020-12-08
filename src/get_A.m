function [A] = get_A(n,flag)
% NOTES from GL: The code is from the test code fd_square.m from the FLAM lib.
% It looks more efficient... but I haven't tried to optimize it.

% initialize
N = (n - 1)^2;  % total number of grid points
switch flag
    case "orig_laplace"
        % set up sparse matrix
        idx = zeros(n+1,n+1);  % index mapping to each point, including "ghost" points
        idx(2:n,2:n) = reshape(1:N,n-1,n-1);
        mid = 2:n;    % "middle" indices -- interaction with self
        lft = 1:n-1;  % "left"   indices -- interaction with one below
        rgt = 3:n+1;  % "right"  indices -- interaction with one above
        I = idx(mid,mid); e = ones(size(I));
        % interactions with ...
        Jl = idx(lft,mid); Sl = -e;                    % ... left
        Jr = idx(rgt,mid); Sr = -e;                    % ... right
        Ju = idx(mid,lft); Su = -e;                    % ... up
        Jd = idx(mid,rgt); Sd = -e;                    % ... down
        Jm = idx(mid,mid); Sm = -(Sl + Sr + Su + Sd);  % ... middle (self)
        % combine all interactions
        I = [ I(:);  I(:);  I(:);  I(:);  I(:)];
        J = [Jl(:); Jr(:); Ju(:); Jd(:); Jm(:)];
        S = [Sl(:); Sr(:); Su(:); Sd(:); Sm(:)];
        % remove ghost interactions
        idx = find(J > 0); I = I(idx); J = J(idx); S = S(idx);
        A = sparse(I,J,S,N,N);
        % clear idx Jl Sl Jr Sr Ju Su Jd Sd Jm Sm I J S
    case "orig_helmholtz"
        error("Not implemented yet!") %TODO
end