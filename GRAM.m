%Implementacao de GramSchimidt 
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


%Implementacao da fatoracao QR de GRAMSCHIMIDT Modificado

function [Q,R] =  mgs(X)
    % Modified Gram-Schmidt.  [Q,R] = mgs(X);
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
    [n,p] = size(X);
    Q = zeros(n,p);
    R = zeros(p,p);
    for k = 1:p
        Q(:,k) = X(:,k);
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(Q(:,k))';
        Q(:,k) = Q(:,k)/R(k,k);
    end
end

[Q,R] = mgs(MATRIZ);


%Resolucao do problema de quadrados minimos com QR
COEF = R \ (Q'*b')



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