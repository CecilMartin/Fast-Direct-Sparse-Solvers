% This script contains a number of functions for testing matrix algebra
% in the "S-matrix" format.
%
% All algorithms are based on recursion on 2x2 block matrices.
%
% NOTE: These codes are intended ONLY to illustrate algorithms; they are
% *exceedingly* inefficient. Both the implementation, and the underlying
% algorithms are highly sub-optimal.
%
% There are three types of matrices:
%   Type 1: Dense with no structure.
%   Type 2: Low-rank.
%   Type 3: 2x2 block matrix (the four block can be any of type 1,2,3).
%
% For a matrix "A", its structured representation "cA" is a "cell"
% structure with entries:
%   cA{1} = type of matrix (1, 2, or 3)
%   cA{2} = size(A,1) - i.e. the number of rows of A
%   cA{3} = size(A,2) - i.e. the number of columns of A
%
% The remaining cells of cA depend on the type of matrix
%
% === Type 1
%   cA{4} = A
%
% === Type 2    A = U * V'
%   cA{4} = U
%   cA{5} = V
%
% === Type 3
%   cA{4} = A11
%   cA{5} = A12
%   cA{6} = A21
%   cA{7} = A22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_S_matrix

rng('default')
rng(0)

LOCAL_verify_add
LOCAL_verify_multiplication
LOCAL_verify_invert
LOCAL_intops

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_verify_add

m = 0.3;
n = 0.6;

%%% Verify 1x1
A  = randn(5,3);
B  = randn(5,3);
cC = LOCAL_add(m,LOCAL_c1(A),n,LOCAL_c1(B));
fprintf(1,'Add 1x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 1x2
A  = randn(5,4);
U  = randn(5,2);
V  = randn(4,2);
B  = U*V';
cC = LOCAL_add(m,LOCAL_c1(A),n,LOCAL_c2(U,V));
fprintf(1,'Add 1x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 1x3
A  = randn(7);
B  = randn(7);
cC = LOCAL_add(m,LOCAL_c1(A),n,LOCAL_c3(B,4));
fprintf(1,'Add 1x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 2x1
U  = randn(7,3);
V  = randn(5,3);
A  = U*V';
B  = randn(7,5);
cC = LOCAL_add(m,LOCAL_c2(U,V),n,LOCAL_c1(B));
fprintf(1,'Add 2x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 2x2
UA = randn(7,3);
VA = randn(5,3);
A  = UA*VA';
UB = randn(7,2);
VB = randn(5,2);
B  = UB*VB';
cC = LOCAL_add(m,LOCAL_c2(UA,VA),n,LOCAL_c2(UB,VB));
fprintf(1,'Add 2x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 2x3
U  = randn(7,3);
V  = randn(7,3);
A  = U*V';
B  = randn(7);
cC = LOCAL_add(m,LOCAL_c2(U,V),n,LOCAL_c3(B,5));
fprintf(1,'Add 2x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 3x1
A  = randn(7);
B  = randn(7);
cC = LOCAL_add(m,LOCAL_c3(A,3),n,LOCAL_c1(B));
fprintf(1,'Add 3x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 3x2
A  = randn(7);
U  = randn(7,2);
V  = randn(7,2);
B  = U*V';
cC = LOCAL_add(m,LOCAL_c3(A,3),n,LOCAL_c2(U,V));
fprintf(1,'Add 3x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

%%% Verify 3x3
A  = randn(7);
B  = randn(7);
cC = LOCAL_add(m,LOCAL_c3(A,3),n,LOCAL_c3(B,3));
fprintf(1,'Add 3x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A - n*B,Inf))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_verify_multiplication

m = 7;

%%% Verify 1x1
A  = randn(5,3);
B  = randn(3,2);
cC = LOCAL_multiply(m,LOCAL_c1(A),LOCAL_c1(B));
fprintf(1,'Multiply 1x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 1x2
A  = randn(5,5);
cA = LOCAL_c1(A);
U  = randn(5,2);
V  = randn(4,2);
B  = U*V';
cB = LOCAL_c2(U,V);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 1x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 1x3
A  = randn(5,7);
cA = LOCAL_c1(A);
B  = randn(7);
cB = LOCAL_c3(B,3);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 1x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 2x1
U  = randn(7,3);
V  = randn(5,3);
A  = U*V';
cA = LOCAL_c2(U,V);
B  = randn(5,2);
cB = LOCAL_c1(B);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 2x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 2x2 (note that there are two versions depending on which k is smaller)
U  = randn(7,3);
V  = randn(5,3);
A  = U*V';
cA = LOCAL_c2(U,V);
UU = randn(5,2);
VV = randn(8,2);
B  = UU*VV';
cB = LOCAL_c2(UU,VV);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 2x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))
U  = randn(7,3);
V  = randn(5,3);
A  = U*V';
cA = LOCAL_c2(U,V);
UU = randn(5,4);
VV = randn(8,4);
B  = UU*VV';
cB = LOCAL_c2(UU,VV);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 2x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 2x3
U  = randn(5,3);
V  = randn(7,3);
A  = U*V';
cA = LOCAL_c2(U,V);
B  = randn(7);
cB = LOCAL_c3(B,4);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 2x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 3x1
A  = randn(7);
cA = LOCAL_c3(A,4);
B  = randn(7,3);
cB = LOCAL_c1(B);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 3x1: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 3x2
A  = randn(7);
cA = LOCAL_c3(A,4);
U  = randn(7,2);
V  = randn(5,2);
B  = U*V';
cB = LOCAL_c2(U,V);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 3x2: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

%%% Verify 3x3
A  = randn(7);
cA = LOCAL_c3(A,4);
B  = randn(7);
cB = LOCAL_c3(B,4);
cC = LOCAL_multiply(m,cA,cB);
fprintf(1,'Multiply 3x3: %12.5e\n',norm(LOCAL_uncompress(cC) - m*A*B,Inf))

return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_verify_invert

acc = 1e-10;

%%% Test just the pure formalism (inversion of 1x1 block matrix).
A  = 7*diag(ones(1,7)) + randn(7);
cB = LOCAL_invert(LOCAL_c1(A),acc);
fprintf(1,'Error in inversion = %12.5e\n',norm(inv(A) - LOCAL_uncompress(cB),Inf))

%%% Test inversion of a 2x2 block matrix.
A  = 7*diag(ones(1,7)) + randn(7);
cB = LOCAL_invert(LOCAL_c3(A,3),acc);
fprintf(1,'Error in inversion = %12.5e\n',norm(inv(A) - LOCAL_uncompress(cB),Inf))

%%% Next attempt inversion of a 4x4 block matrix, where the diagonal blocks
%%% are themselves S-matrices.
A     = 11*diag(ones(1,11)) + randn(11);
cA    = LOCAL_c3(A,5);
cA{4} = LOCAL_c3(LOCAL_uncompress(cA{4}),2);
cA{7} = LOCAL_c3(LOCAL_uncompress(cA{7}),2);
cB    = LOCAL_invert(cA,acc);
fprintf(1,'Error in inversion = %12.5e\n',norm(inv(A) - LOCAL_uncompress(cB),Inf))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is a test program for compression of a matrix emulating the
% behavior of an integral operator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_intops

ntot = 800;
acc  = 1e-10;
nmax = 50;

xx   = sort(rand(1,ntot));
dd   = abs(xx'*ones(1,ntot) - ones(ntot,1)*xx) + eye(ntot);
A    = log(dd);
A    = A + (1+norm(A))*eye(ntot);

cA = LOCAL_compress_brute(A,nmax,acc);
fprintf(1,'Error in compression = %12.5e\n',norm(A - LOCAL_uncompress(cA),Inf))
figure(1)
LOCAL_draw_c3(cA)
title('Compressed version of matrix A')

cB = LOCAL_invert(cA,acc);
fprintf(1,'Error in inversion = %12.5e\n',norm(inv(A) - LOCAL_uncompress(cB),Inf))
figure(2)
LOCAL_draw_c3(cB)
title('Compressed version of matrix B = inv(A)')

%xx   = sort(rand(1,ntot));
%dd   = abs(xx'*ones(1,ntot) - ones(ntot,1)*xx) + eye(ntot);
%B    = log(dd);
%B    = B + (1+norm(B))*eye(ntot);
%cB   = LOCAL_compress_brute(B,nmax,acc);
%cC   = LOCAL_multiply(1,cA,cB);
%fprintf(1,'Error in multiplication = %12.5e\n',norm(A*B - LOCAL_uncompress(cC),Inf))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_draw_c3(cA)

ntot = cA{2};
hold off
plot([0,ntot,ntot,0,0],[0,0,ntot,ntot,0],'k')
hold on
LOCAL_draw_c3_rec(cA,0);
hold off
axis equal
axis ij

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LOCAL_draw_c3_rec(cA,noffset)

n = cA{2};
if (cA{1} == 1)
  plot(noffset + [0.5,n-0.5,n-0.5,0.5,0.5],...
       noffset + [0.5,0.5,n-0.5,n-0.5,0.5],'r')
  text(noffset + 0.5*n, noffset + 0.5*n, sprintf('%d',n))
elseif (cA{1} == 3)
  n1 = cA{4}{2};
  n2 = cA{7}{2};
  %%% Draw the A12 block
  plot(noffset + n1 + [0.5,n2-0.5,n2-0.5,0.5,0.5],...
       noffset +      [0.5,0.5,n1-0.5,n1-0.5,0.5],'g')
  cA12 = cA{5};
  k    = size(cA12{4},2);
  text(noffset + n1 + 0.5*n2, noffset + 0.5*n1, sprintf('%d',k))
  %%% Draw the A12 block
  plot(noffset +      [0.5,n1-0.5,n1-0.5,0.5,0.5],...
       noffset + n2 + [0.5,0.5,n2-0.5,n2-0.5,0.5],'g')
  cA21 = cA{6};
  k    = size(cA21{4},2);
  text(noffset + 0.5*n1, noffset + n1 + 0.5*n2, sprintf('%d',k))
  %%% Draw the two diagonal blocks
  LOCAL_draw_c3_rec(cA{4},noffset)
  LOCAL_draw_c3_rec(cA{7},noffset + cA{4}{2})
else
  fprintf(1,'ERROR: Type 2 matrix on diagonal!\n')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cB = LOCAL_compress_brute(A,nmax,acc)

n = size(A,1);

if (n <= nmax)
  cB = LOCAL_c1(A);
else
  n1 = round(0.5*n + 0.1*(1-2*rand(1)));
  n2 = n - n1;
  I1 = 1:n1;
  I2 = n1 + (1:n2);
  cB = cell(1,7);
  cB{1} = 3;
  cB{2} = n;
  cB{3} = n;
  cB{4} = LOCAL_compress_brute(A(I1,I1),nmax,acc);
  cB{7} = LOCAL_compress_brute(A(I2,I2),nmax,acc);
  %%% Compress A12:
  [U,D,V] = svd(A(I1,I2));
  ss = diag(D);
  k  = sum(diag(D) > acc);
  if (k < 0.25*n)
    VV    = (ones(n2,1)*ss(1:k)').*V(:,1:k);
    cB{5} = LOCAL_c2(U(:,1:k),VV);
  else
    cB{5} = LOCAL_c1(A(I1,I2));
  end
  %%% Compress A21:
  [U,D,V] = svd(A(I2,I1));
  ss = diag(D);
  k  = sum(diag(D) > acc);
  if (k < 0.25*n)
    VV    = (ones(n1,1)*ss(1:k)').*V(:,1:k);
    cB{6} = LOCAL_c2(U(:,1:k),VV);
  else
    cB{6} = LOCAL_c1(A(I2,I1));
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cC = LOCAL_add(m,cA,n,cB)

if (cA{1} == 1)
  if (cB{1} == 1)
    cC = LOCAL_c1(m*cA{4} + n*cB{4});
  elseif (cB{1} == 2)
    cC = LOCAL_c1(m*cA{4} + n*cB{4}*cB{5}');
  elseif (cB{1} == 3)
    cC = LOCAL_c1(m*cA{4} + n*LOCAL_uncompress(cB));
  end
elseif (cA{1} == 2)
  if (cB{1} == 1)
    cC = LOCAL_c1(m*cA{4}*cA{5}' + n*cB{4});
  elseif (cB{1} == 2)
    cC = LOCAL_recompress_c2(LOCAL_c2([m*cA{4},n*cB{4}],[cA{5},cB{5}]),1e-10);
  elseif (cB{1} == 3)
    U     = cA{4};
    V     = cA{5};
    n1    = cB{4}{2};
    n2    = cB{7}{2};
    I1    = 1:n1;
    I2    = n1 + (1:n2);
    cC    = cell(1,7);
    cC{1} = 3;
    cC{2} = cA{2};
    cC{3} = cA{3};
    cA11  = LOCAL_c2(U(I1,:),V(I1,:));
    cA12  = LOCAL_c2(U(I1,:),V(I2,:));
    cA21  = LOCAL_c2(U(I2,:),V(I1,:));
    cA22  = LOCAL_c2(U(I2,:),V(I2,:));
    cC{4} = LOCAL_add(m,cA11,n,cB{4});
    cC{5} = LOCAL_add(m,cA12,n,cB{5});
    cC{6} = LOCAL_add(m,cA21,n,cB{6});
    cC{7} = LOCAL_add(m,cA22,n,cB{7});
  end
elseif (cA{1} == 3)
  if (cB{1} == 1)
    cC = LOCAL_c1(m*LOCAL_uncompress(cA) + n*cB{4});
  elseif (cB{1} == 2)
    U     = cB{4};
    V     = cB{5};
    n1    = cA{4}{2};
    n2    = cA{7}{2};
    I1    = 1:n1;
    I2    = n1 + (1:n2);
    cC    = cell(1,7);
    cC{1} = 3;
    cC{2} = cA{2};
    cC{3} = cA{3};
    cB11  = LOCAL_c2(U(I1,:),V(I1,:));
    cB12  = LOCAL_c2(U(I1,:),V(I2,:));
    cB21  = LOCAL_c2(U(I2,:),V(I1,:));
    cB22  = LOCAL_c2(U(I2,:),V(I2,:));
    cC{4} = LOCAL_add(m,cA{4},n,cB11);
    cC{5} = LOCAL_add(m,cA{5},n,cB12);
    cC{6} = LOCAL_add(m,cA{6},n,cB21);
    cC{7} = LOCAL_add(m,cA{7},n,cB22);
  elseif (cB{1} == 3)
    cC    = cell(1,7);
    cC{1} = 3;
    cC{2} = cA{2};
    cC{3} = cA{3};
    cC{4} = LOCAL_add(m,cA{4},n,cB{4});
    cC{5} = LOCAL_add(m,cA{5},n,cB{5});
    cC{6} = LOCAL_add(m,cA{6},n,cB{6});
    cC{7} = LOCAL_add(m,cA{7},n,cB{7});
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cC = LOCAL_multiply(c,cA,cB)

if (cA{1} == 1)
  if (cB{1} == 1)
    cC = LOCAL_c1(c*cA{4}*cB{4});
  elseif (cB{1} == 2)
    cC = LOCAL_c2(c*cA{4}*cB{4},cB{5});
  elseif (cB{1} == 3)
    A   = cA{4};    % A is a plain matrix.
    n1  = cB{4}{2}; % Extract the nr of rows of the B11-block.
    cA1 = LOCAL_c1(A(:,1:n1));
    cA2 = LOCAL_c1(A(:,(n1+1):end));
    C1  = LOCAL_uncompress(LOCAL_multiply(c,cA1,cB{4})) + ...
          LOCAL_uncompress(LOCAL_multiply(c,cA2,cB{6}));
    C2  = LOCAL_uncompress(LOCAL_multiply(c,cA1,cB{5})) + ...
          LOCAL_uncompress(LOCAL_multiply(c,cA2,cB{7}));
    cC  = LOCAL_c1([C1,C2]);
  end
elseif (cA{1} == 2)
  if (cB{1} == 1)
    cC = LOCAL_c2(c*cA{4},cB{4}'*cA{5});
  elseif (cB{1} == 2)
    kA = size(cA{4},2);
    kB = size(cB{4},2);
    if (kA <= kB)
      cC = LOCAL_c2(c*cA{4},cB{5}*(cB{4}'*cA{5}));
    else
      cC = LOCAL_c2(c*cA{4}*(cA{5}'*cB{4}),cB{5});
    end
  elseif (cB{1} == 3)
    cVt        = LOCAL_c1(cA{5}');
    Vt_times_B = LOCAL_uncompress(LOCAL_multiply(c,cVt,cB));
    cC         = LOCAL_c2(cA{4},Vt_times_B');
  end
elseif (cA{1} == 3)
  if (cB{1} == 1)
    B   = cB{4};     % B is a plain matrix.
    n1  = cA{4}{3};  % Extract the nr of columns of the A11-block
    cB1 = LOCAL_c1(B(1:n1,:));          % Split B = [B1]
    cB2 = LOCAL_c1(B((n1+1):end,:));    %           [B2]
    C1  = LOCAL_uncompress(LOCAL_multiply(c,cA{4},cB1)) + ...
          LOCAL_uncompress(LOCAL_multiply(c,cA{5},cB2));
    C2  = LOCAL_uncompress(LOCAL_multiply(c,cA{6},cB1)) + ...
          LOCAL_uncompress(LOCAL_multiply(c,cA{7},cB2));
    cC  = LOCAL_c1([C1;C2]);
  elseif (cB{1} == 2)
    cU = LOCAL_c1(cB{4});
    A_times_U = LOCAL_uncompress(LOCAL_multiply(c,cA,cU));
    cC = LOCAL_c2(A_times_U,cB{5});
  elseif (cB{1} == 3)
    cC    = cell(1,7);
    cC{1} = 3;
    cC{2} = cA{4}{2} + cA{6}{2};
    cC{3} = cB{4}{3} + cB{5}{3};
    cC{4} = LOCAL_add(1,LOCAL_multiply(c,cA{4},cB{4}),1,LOCAL_multiply(c,cA{5},cB{6}));
    cC{5} = LOCAL_add(1,LOCAL_multiply(c,cA{4},cB{5}),1,LOCAL_multiply(c,cA{5},cB{7}));
    cC{6} = LOCAL_add(1,LOCAL_multiply(c,cA{6},cB{4}),1,LOCAL_multiply(c,cA{7},cB{6}));
    cC{7} = LOCAL_add(1,LOCAL_multiply(c,cA{6},cB{5}),1,LOCAL_multiply(c,cA{7},cB{7}));
    fprintf(1,'WARNING: After multiplying type 3 by type 3, recompression should be done for any type-2 matrix.\n')
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function exploits a formula for the inverse of a 2x2 block matrix;
%
%    inv([A11 A12]) = [B11 B12]
%        [A21 A22]    [B21 B22]
%
% where
%
%    B11 = inv(A11 - A12 * inv(A22) * A21)
%    B12 = -B11 * A12 * inv(A22)
%    B21 = -inv(A22) * A21 * B11
%    B22 = inv(A22) + inv(A22) * A21 * B11 * A12 * inv(A22)
%
% The point is that the formulas require two inverses of matrices half the size.
% Since many groups recur, we introduce intermediate variables as follows:
%
%    X22 = inv(A22)
%    X12 = A12 * inv(A22) = A12 * X22
%    X21 = inv(A22) * A21 = X22 * A21
% 
% Then
%
%    B11 = inv(A11 - X12 * A21)
%    B12 = -B11 * X12
%    B21 = -X21 * B11
%    B22 = X22 - B21 * X12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cB = LOCAL_invert(cA,acc)

if (cA{1} == 1)
  cB = LOCAL_c1(inv(cA{4}));
elseif (cA{1} == 2)
  fprintf(1,'ERROR: You are trying to invert a low-rank matrix.\n')
  keyboard
elseif (cA{1} == 3)
  cX22  = LOCAL_invert(cA{7},acc);        % X22  = inv(A22)
  cX12  = LOCAL_multiply(1,cA{5},cX22);   % X12  = A12 * X22
  cX21  = LOCAL_multiply(1,cX22,cA{6});   % X21  = X22 * A21
  cTMP  = LOCAL_multiply(1,cX12,cA{6});   % TMP  = X12 * A21
  cTMP2 = LOCAL_add(1,cA{4},-1,cTMP);     % TMP2 = A11 - TMP
  cB11  = LOCAL_invert(cTMP2,acc);        % B11  = inv(TMP2)
  cB12  = LOCAL_multiply(-1,cB11,cX12);   % B12  = B11 * X12
  cB21  = LOCAL_multiply(-1,cX21,cB11);   % B21  = X21 * B11
  cTMP  = LOCAL_multiply(1,cB21,cX12);    % TMP  = B21 * X22
  cB22  = LOCAL_add(1,cX22,-1,cTMP);      % B22  = X22 - TMP
  cB    = cell(1,7);
  cB{1} = 3;
  cB{2} = cA{2};
  cB{3} = cA{3};
  cB{4} = cB11;
  if (cB12{1} == 2)
    cB{5} = LOCAL_recompress_c2(cB12,acc);
  else
    cB{5} = cB12;
  end
  if (cB21{1} == 2)
    cB{6} = LOCAL_recompress_c2(cB21,acc);
  else
    cB{6} = cB21;
  end
  cB{7} = cB22;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cB = LOCAL_recompress_c2(cA,acc)

if ~(cA{1} == 2)
  fprintf(1,'ERROR: The input is not of type 2.\n')
  keyboard
end

%%% Extract the factors in the low-rank factorization.
U = cA{4};
V = cA{5};

%%% Orthonormalize the columns of U:
[Q,R] = qr(U,0);

%%% Compress the matrix R*V' down to rank k and form the new factors.
[VV,DD,UU] = svd(V*R',0);
k          = sum(diag(DD) >= acc);
UNEW       = Q*UU(:,1:k);
VNEW       = VV(:,1:k)*DD(1:k,1:k);

%%% Store the new factors away.
cB = LOCAL_c2(UNEW,VNEW);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_uncompress(cA)

if (cA{1} == 1)
  A = cA{4};
elseif (cA{1} == 2)
  A = cA{4}*cA{5}';
elseif(cA{1} == 3)
  A11 = LOCAL_uncompress(cA{4});
  A12 = LOCAL_uncompress(cA{5});
  A21 = LOCAL_uncompress(cA{6});
  A22 = LOCAL_uncompress(cA{7});
  A   = [A11,A12;...
         A21,A22];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cA = LOCAL_c1(A)
cA    = cell(1,4);
cA{1} = 1;
cA{2} = size(A,1);
cA{3} = size(A,2);
cA{4} = A;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cA = LOCAL_c2(U,V)
cA    = cell(1,5);
cA{1} = 2;
cA{2} = size(U,1);
cA{3} = size(V,1);
cA{4} = U;
cA{5} = V;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cA = LOCAL_c3(A,n1)

ntot  = size(A,1);
I1    = 1:n1;
I2    = (n1+1):ntot;
cA    = cell(1,7);
cA{1} = 3;
cA{2} = ntot;
cA{3} = ntot;
cA{4} = LOCAL_c1(A(I1,I1));
cA{5} = LOCAL_c1(A(I1,I2));
cA{6} = LOCAL_c1(A(I2,I1));
cA{7} = LOCAL_c1(A(I2,I2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

