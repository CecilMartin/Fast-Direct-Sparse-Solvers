% main
n1 = 16;
n2 = 16;
N = n1*n2;
n=[n1,n2];
D = 2; %2D
flag = "orig_laplace";
% orig_laplace, standard 5pt laplace
% orig_helmholtz, standart 5pt laplace + centered difference for 1st order
% term
BV_flag = "Dirichlet"; % "Neumann"
ordering_flag = "Nested"; % "Sweeping"

A = get_A(n,flag,BV_flag);
ind = ordering(A,X,ordering_flag);
L,U = my_solve(A,X,ind,.......);



