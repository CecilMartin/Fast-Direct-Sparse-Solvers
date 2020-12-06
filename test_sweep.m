% NOTES from GL: I've done the first part of the sweeping alogrithm (the
% computation of Schur complement part). One of the current problems is
% that now the size n (number of points for one column) have to be a square
% number, (like to say, 64 is okay but 128 is not). That is mainly due to
% I'm a little confused on how to use the function 'hifie2', which is used
% to compute the 'efficient' factorization of the given matrix. 

n = 4^4;  % Number of points for one column, total points would be n^2
occ = 4;  % Parameter for the factorization, looks like the size of matrix on the lead node
rank_or_tol = 1e-9; % Tolerance for the rank approximation (epsilon)

[A] = get_A_orig_test(n+1); % Get the sparse stiffness matrix A

tic;
I1 = 1:n;
% x is an important input for hifie2, but I'm not sure the meaning of it
[x1,x2] = ndgrid((1:sqrt(n))/sqrt(n)); x = [x1(:) x2(:)]';
% Efficient factorization and record it in the cell
F_{1} = hifie2(A(I1, I1),x,occ,rank_or_tol,[], []);
for i = 2:n
    Ii = (i-1)*n + I1;
    % Compute the Schur complement (I didn't do the inverse, the Xk as the 
    % algorithm in the book, but merely the Schur complement Sk)
    S = A(Ii, Ii) - A(Ii, Ii-n) * hifie_sv(F_{i-1}, A(Ii-n, Ii));
    % Efficient factorization and record it in the cell
    F_{i} = hifie2(S,x,occ,rank_or_tol,[], []);
end
t = toc;
fprintf('hifie2 time: %10.4e (s) \n',t);