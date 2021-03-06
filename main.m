% NOTES from GL: I've done the first part of the sweeping alogrithm (the
% computation of Schur complement part). One of the current problems is
% that now the size n (number of points for one column) have to be a square
% number, (like to say, 64 is okay but 128 is not). That is mainly due to
% I'm a little confused on how to use the function 'hifie2', which is used
% to compute the 'efficient' factorization of the given matrix. 

% The real solution is u = sin(pi*x) * cos(pi*y)
% Equation: - Delta u = f = 2*pi^2 u = 2*pi^2 * sin(pi*x) * cos(pi*y)
clear;close all;clc;

% PATH
addpath(genpath(pwd));
% FLAM_path = '/home/cecil/Project/FLAM';
FLAM_path = '../../Final project/Code/FLAM-master';
if ~(exist('hifie2'))
    current_path = pwd;
    cd(FLAM_path);
    startup;
    cd(current_path);
end

D = 2; %2D
% flag = "orig_laplace";
% orig_laplace, standard 5pt laplace
flag = "orig_helmholtz";
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"


n = 64 ;  % Number of points for one column, total points would be n^2
occ = 8;  % Parameter for the factorization, looks like the size of matrix on the lead node
rank_or_tol = 1e-9; % Tolerance for the rank approximation (epsilon)

kappa = 10;
f_fun_test = @(x) (2*pi^2 - kappa^2) * sin(pi*x(:, 1)) .* cos(pi*x(:, 2));

[A] = get_A(n+1,flag, kappa); % Get the sparse stiffness matrix A
% Get the load matrix
% f = randn(n^2,1); 

f = get_f(n+1,@bv_fun_test,f_fun_test,flag,kappa);

tic;

u = sweeping(n,A,f,occ,rank_or_tol);

t = toc;
fprintf('hifie2 time: %10.4e (s) \n',t);

[x1, x2] = ndgrid((1:n)/(n+1)); x = [x1(:) x2(:)];
u_rel = u_fun_test(x);

fprintf('Error of the solution: %10.4e \n', norm(u - u_rel) / norm(u_rel));