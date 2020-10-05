%Implementacao de Householder
%Algoritmo de fatoração QR foi pego numa publicação de
%  Cleve Moler, July 25, 2016
%Com referencias para G. W. (Pete) Stewart
%referencias para publicação completa na bibliografia do projeto

%Matriz de vader
p = fliplr(vander([0, 1/20, 2/20, 3/20, 4/20, 5/20, 6/20, 7/20, 8/20, 9/20, 10/20, 11/20, 12/20, 13/20, 14/20, 15/20, 16/20, 17/20, 18/20, 19/20, 1]));
MATRIZ = p(:,[1:11]);

%F(x)

b = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]*0.05;
b = exp(sin(6*b));


%Implementacao da fatoracao QR de HOUSEHOLDER
% This program does not actually compute the QR orthogonalization, but rather computes 
%R and a matrix U containing vectors that generate the Householder reflectors 
%whose product is Q.
function [u,nu] = housegen(x)
    % [u,nu] = housegen(x)
    % Generate Householder reflection.
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    % [u,nu] = housegen(x).
    % H = I - uu' with Hx = -+ nu e_1
    %    returns nu = norm(x).
    u = x;
    nu = norm(x);
    if nu == 0
        u(1) = sqrt(2);
        return
    end
    u = x/nu;
    if u(1) >= 0
        u(1) = u(1) + 1;
        nu = -nu;
    else
        u(1) = u(1) - 1;
    end
    u = u/sqrt(abs(u(1)));
end

function [U,R] = hqrd(X)  
    % Householder triangularization.  [U,R] = hqrd(X);
    % Generators of Householder reflections stored in U.
    % H_k = I - U(:,k)*U(:,k)'.
    % prod(H_m ... H_1)X = [R; 0]
    % where m = min(size(X))
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    U = zeros(size(X));
    m = min(n,p);
    R = zeros(m,m);
    for k = 1:min(n,p)
        [U(k:n,k),R(k,k)] = housegen(X(k:n,k));
        v = U(k:n,k)'*X(k:n,k+1:p);
        X(k:n,k+1:p) = X(k:n,k+1:p) - U(k:n,k)*v;
        R(k,k+1:p) = X(k,k+1:p);
    end
end


[U,R] =hqrd(MATRIZ);

function Z = house_apply(U,X)
    % Apply Householder reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q*X without actually computing Q.
    H = @(u,x) x - u*(u'*x);
    Z = X;
    [~,n] = size(U);
    for j = n:-1:1
        Z = H(U(:,j),Z);
    end
end

function Z = house_apply_transpose(U,X)
    % Apply Householder transposed reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q'*X without actually computing Q'.
    H = @(u,x) x - u*(u'*x);
    Z = X;
    [~,n] = size(U);
    for j = 1:n
        Z = H(U(:,j),Z);
    end
end

I = eye(size(U));
   Q = house_apply(U,I)
%contrucao da Q 
 % H_k = I - U(:,k)*U(:,k)'.
 %

%Resolucao do problema de quadrados minimos com QR
COEF = R \ (Q'*b')(1:11)

%obtem valores do polinomio aplicado nos x abaixo
X = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]*0.05;

for j = 1:21
  x = X(j);
DIAG  = zeros(11);
for i = 1:11
    DIAG(i,i) = x^(i-1);
end
z(j) = sum(COEF'*DIAG);
end
z'
%calcula diferenca entre valor de f(x) e valor obtido por polinomio
D = b'-z'
%calcula numeros de condicao
CQ = cond(Q,2)
CR = cond(R,2)

COEF


