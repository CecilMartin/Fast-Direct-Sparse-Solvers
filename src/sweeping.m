function u = sweeping(n,A,f,occ,rank_or_tol)
I1 = 1:n;
% x is an important input for hifie2, but I'm not sure the meaning of it
sn = round(sqrt(n));
[x1,x2] = ndgrid((1:sn)/sn); x = [x1(:) x2(:)]';
opts = struct('Tmax',2,'skip',1,'verb',0);
% Efficient factorization and record it in the cell
F_all{1} = hifie2(A(I1, I1),x,occ,rank_or_tol,[], opts);
f_tilde = zeros(n);
f_tilde(:, 1) = f(I1);
for i = 2:n
    Ii = (i-1)*n + I1;
    % Compute the Schur complement (I didn't do the inverse, the Xk as the 
    % algorithm in the book, but merely the Schur complement Sk)
    S = A(Ii, Ii) - A(Ii, Ii-n) * hifie_sv(F_all{i-1}, A(Ii-n, Ii)); % Sparse-dense multi

    % Efficient factorization and record it in the cell
    F_all{i} = hifie2(S,x,occ,rank_or_tol,[], opts);
    f_tilde(:, i) = f(Ii) - A(Ii, Ii-n)*hifie_sv(F_all{i-1}, f_tilde(:, i-1));
    
end



u = zeros(size(f));
u(Ii) = hifie_sv(F_all{n}, f_tilde(:, n));
for i = (n-1):-1:1
    Ii = (i-1)*n+I1;
    u(Ii) = hifie_sv(F_all{i}, f_tilde(:, i)-A(Ii,Ii+n)*u(Ii+n));
end
end