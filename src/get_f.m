function f = get_f(n,bv_fun,f_fun,flag)
% load vector for Laplace
% Dirichlet Boundary Value
% The region is set as [0, 1] * [0, 1]
% With (n-1)*(n-1) interior points
% f = -ones((n-1)^2,1);
h = 1/n;

[x1, x2] = ndgrid((1:(n-1))*h); x = [x1(:) x2(:)];
f = h^2 * f_fun(x); %TODO

switch flag
    case "orig_laplace"
        Il = 1:n-1; x_tmp = [x(Il,1),x(Il,2)-h]; 
        f(Il) = f(Il) + bv_fun(x_tmp);  % The left boundary, {0}*[0, 1]

        Ir = (n-1)*(n-2) + (1:n-1); x_tmp = [x(Ir,1),x(Ir,2)+h]; 
        f(Ir) = f(Ir) + bv_fun(x_tmp);  % The right boundary, {1}*[0, 1]

        Iu = 1:(n-1):(n-1)^2-(n-2); x_tmp = [x(Iu,1)-h,x(Iu,2)]; 
        f(Iu) = f(Iu) + bv_fun(x_tmp);  % The upper boundary, [0, 1]*{1}

        Id = (n-1):(n-1):(n-1)^2; x_tmp = [x(Id,1)+h,x(Id,2)]; 
        f(Id) = f(Id) + bv_fun(x_tmp);  % The down boundary, [0, 1]*{0}
    case "orig_helmholtz"
        error("Not implemented yet") %TODO
end