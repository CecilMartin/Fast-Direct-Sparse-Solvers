function u = sweeping_buffered(n,m,b,A,f,occ,rank_or_tol)
assert(n==1+m*(b+1))

% Update the buffered stiffness matrix and source term
A_tilde = sparse((m+1)*n, (m+1)*n);
f_tilde = zeros((m+1)*n, 1);
I1 = 1:n;

for k = 1:(m+1)
    I2km1 = (k-1)*(b+1)*n + I1;
    A_tilde((k-1)*n + I1, (k-1)*n + I1) = A(I2km1, I2km1);
    f_tilde((k-1)*n + I1) = f(I2km1);
end

for k = 1:m
    I2km1 = (k-1)*(b+1)*n + I1;  % Index of I_2k-1 in original A
    I2kp1 = k*(b+1)*n + I1;  % Index of I_2k+1 in original A
    I2k = (k-1)*(b+1)*n + n + (1:(n*b));
    A_tilde((k-1)*n + I1, (k-1)*n + I1) = ...
        A_tilde((k-1)*n + I1, (k-1)*n + I1) - A(I2km1, I2k) * (A(I2k, I2k) \ A(I2k, I2km1));
    A_tilde(k*n + I1, k*n + I1) = ...
        A_tilde(k*n + I1, k*n + I1) - A(I2kp1, I2k) * (A(I2k, I2k) \ A(I2k, I2kp1));
    A_tilde((k-1)*n + I1, k*n + I1) = - A(I2km1, I2k) * (A(I2k, I2k) \ A(I2k, I2kp1));
    A_tilde(k*n + I1, (k-1)*n + I1) = - A(I2kp1, I2k) * (A(I2k, I2k) \ A(I2k, I2km1));
    f_tilde((k-1)*n + I1) = ...
        f_tilde((k-1)*n + I1) - A(I2km1, I2k) * (A(I2k, I2k) \ f(I2k));
    f_tilde(k*n + I1) = ...
        f_tilde(k*n + I1) - A(I2kp1, I2k) * (A(I2k, I2k) \ f(I2k));
end


I1 = 1:n;
opts = struct('Tmax',2,'skip',1,'verb',0);
% x is an important input for hifie2, but I'm not sure the meaning of it
sn = round(sqrt(n));
[x1,x2] = ndgrid((1:sn)/sn); x = [x1(:) x2(:)]';
% Efficient factorization and record it in the cell
F_all{1} = hifie2(A_tilde(I1, I1),x,occ,rank_or_tol,[], opts);
for k = 1:m
    Ik = k*n + I1;
    % Compute the Schur complement (I didn't do the inverse, the Xk as the 
    % algorithm in the book, but merely the Schur complement Sk)
    S = A_tilde(Ik, Ik) - A_tilde(Ik, Ik-n) * hifie_sv(F_all{k}, A_tilde(Ik-n, Ik)); % Sparse-dense multi
    % Efficient factorization and record it in the cell
    F_all{k+1} = hifie2(S,x,occ,rank_or_tol,[], opts);   
    f_tilde(Ik) = f_tilde(Ik) - A_tilde(Ik, Ik-n)*hifie_sv(F_all{k}, f_tilde(Ik-n));
end


% Solve for u on survival columns
u = zeros(size(f));
I2mp1 = m*(b+1)*n + I1;
u(I2mp1) = hifie_sv(F_all{m+1}, f_tilde(Ik));
for k = m:-1:1
    Ik = (k-1)*n + I1;
    I2km1 = (k-1)*(b+1)*n + I1;
    u(I2km1) = hifie_sv(F_all{k}, f_tilde(Ik)-A_tilde(Ik,Ik+n)*u(I2km1+(b+1)*n));
end

% Solve for u on buffered columns
for i = 1:m
    I2i = (i-1)*(b+1)*n + n + (1:(n*b));
    I2im1 = (i-1)*(b+1)*n + I1;  % Index of I_2k-1 in original A
    I2ip1 = i*(b+1)*n + I1;  % Index of I_2k+1 in original A
    u(I2i) = A(I2i, I2i) \ (f(I2i) - A(I2i, I2im1)*u(I2im1) - A(I2i, I2ip1)*u(I2ip1));
end
end