%Algoritimo da fatoração de Cholesky foi adaptado do 
%Algoritimo da Profa. Dra. Marli de Freitas Gomes Hernandez
%CESET-UNICAMP
%dados dos experimentos no Si, D em metros e T em segundos
D =[ 0.896,0.896,0.896,0.796,0.796,0.796,0.696, 0.696,0.696,0.596,0.596,0.596,0.496,0.496,0.496,0.396, 0.396,0.396];
T = [2.076,2.084,2.082,1.831,1.827,1.937,1.957,1.958,1.951,1.905,1.910,1.908,1.883,1.884,1.877,1.900,1.901,1.902];
D = D';
T=T';
A=ones(18,2);
C = ones(18,1);
for i = 1:18
  A(i,1) = D(i)*D(i);
  C(i) = T(i)*T(i)*D(i);
end


%eq normais AtAx = Atb 
ATA = A'*A;
ATB = A'* C;


%cholesky

n =2;

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

COEF = r\y;

g = 4*pi*pi/COEF(1)
k = sqrt(COEF(2)*g)
