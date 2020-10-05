%Algoritimo da fatoração de Cholesky foi adaptado do 
%Algoritimo da Profa. Dra. Marli de Freitas Gomes Hernandez
%CESET-UNICAMP

%Matriz de vader
p = fliplr(vander([0, 1/20, 2/20, 3/20, 4/20, 5/20, 6/20, 7/20, 8/20, 9/20, 10/20, 11/20, 12/20, 13/20, 14/20, 15/20, 16/20, 17/20, 18/20, 19/20, 1]));
MATRIZ = p(:,[1:11]);

%F(x)

b = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]*0.05;
b = exp(sin(6*b));
b = b';

%eq normais AtAx = Atb 
ATA = MATRIZ'*MATRIZ;
ATB = MATRIZ'* b;


%cholesky

n = 11;

A = ATA;
% Determinar a Matriz R tal que A=(R�)R, A sim�trica positiva definida


% testar simetria da matriz A
for i=1:n
    for j=i+1:n
        if(A(i,j)~=A(j,i))
            disp('ERRO: matriz A assim�trica')
            c=input('digite CONTROL^C  '); % forma de parar o programa
        end
    end
end
% Determinar G tal que A=G'G
if A(1,1) < 0 
    disp('ERRO: raiz quadrada de numero negativo g(i,i).')
    c=input('digite CONTROL^C  '); % forma de parar o programa
else
   if A(1,1) == 0
        disp('ERRO: divis�o por 0 no c�lculo de g[i,j]/g(i,i).')
        c=input('digite CONTROL^C  '); % forma de parar o programa
   else
        g(1,1)=sqrt(A(1,1));
        for j=2:n
            g(1,j)= A(1,j)/g(1,1);
        end
        for i=2:n
            g(i,i)=A(i,i);
            for k= 1:i-1
                g(i,i) =g(i,i) - g(k,i)^2;
            end
            if g(i,i) < 0 
                disp('ERRO: raiz quadrada de numero negativo g(i,i).')
                c=input('digite CONTROL^C'); % forma de parar o programa
            else
                 if g(i,i) == 0
                     disp('ERRO: divis�o por 0 no c�lculo de g[i,j]/g(i,i).')
                     c=input('digite CONTROL^C'); % forma de parar o programa
                 else
                 g(i,i)=sqrt(g(i,i));
                 for j=i+1:n
                    g(i,j)=A(i,j);
                    for k=1:i-1
                        g(i,j)=g(i,j)-g(k,i)*g(k,j);
                    end
                    g(i,j)=g(i,j)/g(i,i);
                  end
                end
            end
        end
   end
end
% Outra forma de programar o Cholesky
% Determinar R tal que A=R'R
for i=1:n
  for j=i:n
       r(i,j) = A(i,j);
       if i == j 
           if i > 1
               for k=1:i-1
                   r(i,j) = r(i,j) - r(k,i)^2; 
               end
           end
           if r(i,i) < 0 
                disp('ERRO: raiz quadrada de numero negativo g(i,i).')
                c=input('digite CONTROL^C'); % forma de parar o programa
           else
              if r(i,i) == 0
                  disp('ERRO: divis�o por 0 no c�lculo de g[i,j]/g(i,i).')
                  c=input('digite CONTROL^C'); % forma de parar o programa
              else
                  r(i,j)=sqrt(r(i,j));
              end
           end
        else
           if i > 1
               for k=1:i-1
                   r(i,j) = r(i,j) - r(k,i)*r(k,j); 
               end
           end
           r(i,j) = r(i,j)/r(i,i);
       end
  end
end


%ATA = R'R 
%R'Y=AtB
y = r'\ATB;

%RCOEF = Y

COEF = r\y

%obtem valores do polinomio aplicado nos x abaixo (z = valores nos polinomios)
X = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]*0.05;

for j = 1:21
  x = X(j);
DIAG  = zeros(11);
for i = 1:11
    DIAG(i,i) = x^(i-1);
end
z(j) = sum(COEF'*DIAG);
end


%calcula diferenca entre valor de f(x) e valor obtido por polinomio
D = b-z'


 %calcula numero de condicao 
CR = cond(r,inf)
COEF 
