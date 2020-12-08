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

%% Generate Figure 1
% Computing time of Laplacian Eq vs n, different epsilon

D = 2; %2D
flag = "orig_laplace";
% orig_laplace, standard 5pt laplace
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"

tol_lst = [1e-6, 1e-9, 1e-12];
n_lst = round(sqrt(2.^(2:7))).^2;

occ = 8;  % Parameter for the factorization, looks like the size of matrix on the lead node

time_orig = zeros(length(tol_lst), length(n_lst));
for m = 1:length(tol_lst) % Tolerance for the rank approximation (epsilon)
    tol = tol_lst(m);
    lgd{m} = sprintf('tol = %2.1e', tol);
    for k = 1:length(n_lst) % Number of points for one column, total points would be n^2
        n = n_lst(k);
        [A] = get_A(n+1,flag); % Get the sparse stiffness matrix A
        f = get_f(n+1,@bv_fun_test,@f_fun_test,flag); % Get the load matrix
        tic;
        u = sweeping(n,A,f,occ,tol);
        t = toc;

        [x1, x2] = ndgrid((1:n)/(n+1)); x = [x1(:) x2(:)];
        u_rel = u_fun_test(x);
        time_orig(m, k) = t;
        % fprintf('tol = %2.1e, n = %d finished! \n', tol, n);
        fprintf('tol = %2.1e, n = %d, time: %10.4e (s), Error of the solution: %10.4e \n', tol, n, t, norm(u - u_rel) / norm(u_rel));
    end
end

fig1 = figure(1);
loglog(n_lst, time_orig, 's-');
xlabel('n (number of points of one line)');
ylabel('Time');
title('Running time of Laplacian Eq');
legend(lgd, 'Location', 'northwest');
set(gca,'FontSize',16);


%% Generate Figure 2
% Computing time of Helmholtz Eq vs n, different kappa
clear;

D = 2; %2D
flag = "orig_helmholtz";
% orig_laplace, standard 5pt laplace
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"

rank_or_tol = 1e-9;
kappa_lst = [0, 1, 10, 20];
n_lst = round(sqrt(2.^(2:7))).^2;

occ = 8;  % Parameter for the factorization, looks like the size of matrix on the lead node

time_orig = zeros(length(kappa_lst), length(n_lst));
for m = 1:length(kappa_lst) % Tolerance for the rank approximation (epsilon)
    kappa = kappa_lst(m);
    lgd{m} = sprintf('kappa = %d', kappa);
    for k = 1:length(n_lst) % Number of points for one column, total points would be n^2
        n = n_lst(k);
        [A] = get_A(n+1,flag,kappa); % Get the sparse stiffness matrix A
        u_fun_helm = @(x) sin(kappa*x(:, 1)) + cos(kappa*x(:, 2));
        bv_fun_helm = @(x) u_fun_helm(x);
        f = get_f(n+1,bv_fun_helm,@f_fun_test,flag,kappa); % Get the load matrix
        tic;
        u = sweeping(n,A,f,occ,rank_or_tol);
        t = toc;

        [x1, x2] = ndgrid((1:n)/(n+1)); x = [x1(:) x2(:)];
        u_rel = u_fun_helm(x);
        time_orig(m, k) = t;
        % fprintf('kappa = %d, n = %d finished! \n', kappa, n);
        fprintf('kappa = %d, n = %d, time: %10.4e (s), Error of the solution: %10.4e \n', kappa, n, t, norm(u - u_rel) / norm(u_rel));
    end
end
fig2 = figure(2);
loglog(n_lst, time_orig, 's-');
xlabel('n (number of points of one line)');
ylabel('Time');
title('Running time of Laplacian Eq');
legend(lgd, 'Location', 'northwest');
set(gca,'FontSize',16);

%% Generate Figure 3
% Computing time of Laplacian Eq vs different buffer size
clear;

D = 2; %2D
flag = "orig_laplace";
% orig_laplace, standard 5pt laplace
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"

tol_lst = [1e-6, 1e-9, 1e-12];
% n = 4^4; b_lst = [0, 2, 4, 14, 16, 50]; 
n = 13^2; b_lst = [0, 1, 2, 3, 5, 6, 7, 11, 20, 23, 41, 83];

occ = 8;  % Parameter for the factorization, looks like the size of matrix on the lead node
rank_or_tol = 1e-9;

time_buffer = zeros(length(b_lst), 1);
for k = 1:length(b_lst) % Tolerance for the rank approximation (epsilon)
    b = b_lst(k);
    m = round((n-1)/(b+1));
    if abs(1+m*(b+1) - n) > 0.5
        error('m, b not legal \n');
    end
    [A] = get_A(n+1,flag); % Get the sparse stiffness matrix A
    f = get_f(n+1,@bv_fun_test,@f_fun_test,flag); % Get the load matrix
    tic;
    if b == 0
        u = sweeping(n,A,f,occ,rank_or_tol);
    else
        u = sweeping_buffered(n,m,b,A,f,occ,rank_or_tol);
    end
    t = toc;

    [x1, x2] = ndgrid((1:n)/(n+1)); x = [x1(:) x2(:)];
    u_rel = u_fun_test(x);
    time_buffer(k) = t;
    % fprintf('b = %d, m = %d finished! \n', b, m);
    fprintf('b = %d, m = %d, time: %10.4e (s), Error of the solution: %10.4e \n', b, m, t, norm(u - u_rel) / norm(u_rel));
end

fig3 = figure(3);
plot(b_lst, time_buffer, 's-');
xlabel('b (width of buffer)');
ylabel('Time');
title('Running time of Laplacian Eq');
set(gca,'FontSize',16);