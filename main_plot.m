clear;close all;clc;

set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',22)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

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

tol_lst = [1e-3, 1e-8, 1e-13];
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
loglog(n_lst.^2, time_orig, 's-');
xlabel('N');
ylabel('Time(s)');
% title('Running time of Laplacian Eq');
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

rank_or_tol = 1e-12;
kappa_lst = [0, 1, 10, 20];
n_lst = round(sqrt(2.^(2:7))).^2;

occ = 4;  % Parameter for the factorization, looks like the size of matrix on the lead node

time_orig = zeros(length(kappa_lst), length(n_lst));
for m = 1:length(kappa_lst) % Tolerance for the rank approximation (epsilon)
    kappa = kappa_lst(m); 
    lgd{m} = ['$\kappa$ = 2 $\pi$ * ',num2str(kappa)];
    kappa = kappa * (1*pi);
    for k = 1:length(n_lst) % Number of points for one column, total points would be n^2
        n = n_lst(k);
        [A] = get_A(n+1,flag,kappa); % Get the sparse stiffness matrix A
%         u_fun_helm = @(x) sin(kappa/sqrt(2)*x(:, 1)) .* cos(kappa/sqrt(2)*x(:, 2));
        u_fun_helm = @(x) sin(kappa*x(:, 1)) + cos(kappa*x(:, 2)) ...
            + sin(pi*x(:, 1)) .* exp(pi*x(:, 2));
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
loglog(n_lst.^2, time_orig, 's-');
xlabel('N');
ylabel('Time');
% title('Running time of Laplacian Eq');
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
% title('Running time of Laplacian Eq');
set(gca,'FontSize',16);

%% Generate Figure 4
% Skeletonization factorization ID

D = 2; %2D
flag = "orig_laplace";
% orig_laplace, standard 5pt laplace
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"

tol_lst = [1e-3, 1e-6, 1e-9, 1e-12, 1e-14];
n_lst = 4.^(1:4);

occ = 16;  % Parameter for the factorization, looks like the size of matrix on the lead node
repeat = 8;

time_orig = zeros(length(tol_lst), length(n_lst));
for m = 1:length(tol_lst) % Tolerance for the rank approximation (epsilon)
    tol = tol_lst(m);
    lgd{m} = sprintf('tol = %2.1e', tol);
    for k = 1:length(n_lst) % Number of points for one column, total points would be n^2
        n = n_lst(k);
        [A] = get_A(n+1,flag); % Get the sparse stiffness matrix A
        f = get_f(n+1,@bv_fun_test,@f_fun_test,flag); % Get the load matrix
        I1 = 1:n;
        % x is an important input for hifie2, but I'm not sure the meaning of it
        sn = round(sqrt(n));
        [x1,x2] = ndgrid((1:sn)/sn); x = [x1(:) x2(:)]';
        opts = struct('Tmax',2,'skip',1,'verb',0);
        % Pre-run       
        F_tmp{1} = hifie2(A(I1, I1),x,occ,tol,[], opts);
        for i = 2:4
            Ii = (i-1)*n + I1;
            S = A(Ii, Ii) - A(Ii, Ii-n) * hifie_sv(F_tmp{i-1}, A(Ii-n, Ii)); % Sparse-dense multi

            % Efficient factorization and record it in the cell
            F_tmp{i} = hifie2(S,x,occ,tol,[], opts);
        end
        tic;
        for i = 1:repeat
            tmp = hifie2(A,n+1,occ,tol, [], opts);
        end
        t = toc;

%         [x1, x2] = ndgrid((1:n)/(n+1)); x = [x1(:) x2(:)];
%         u_rel = u_fun_test(x);
        time_orig(m, k) = t/repeat;
        fprintf('tol = %2.1e, n = %d finished! \n', tol, n);
%         fprintf('tol = %2.1e, n = %d, time: %10.4e (s), Error of the solution: %10.4e \n', tol, n, t, norm(u - u_rel) / norm(u_rel));
    end
end

fig1 = figure(4);
loglog(n_lst.^2, time_orig, 's-');
xlabel('N');
ylabel('Time');
% title('Running time of Laplacian Eq');
legend(lgd, 'Location', 'northwest');
set(gca,'FontSize',16);

