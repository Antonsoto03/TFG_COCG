
%A = mmread('young2c.mtx'); %USANDO LA FUNCIÓN mmread.m
%A = [  3,      1+2i,    0.5-1i,  -2i ;
 %1+2i,   2,       4i,      1-1i;
  %     0.5-1i, 4i,     -1,       3+2i;
   %   -2i,     1-1i,    3+2i,    0    ]; 
   %introduciendo la matriz manualmente

%%%%%%%%%%%%%%%%%%%%%
n = 1200; %Número de columnas de la matriz
rng(0)
A = rnd_csPD(n,0.4); %matriz generada aleatoriamente
b = randn(n,1)+1i*randn(n,1); %termino independiente generado aleatoriamente
%b= ones(n,1)+1i*ones(n,1);
x0 = zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%cond(A); % la funcion cond no funciona con matrices sparse
%b = ones(841,1) ; % Vector termino independiente
%b = [ 1; 2+ 1i;-1i;3-2i ]; %Vector termino independiente introducido
%manualmente
tol = 1e-8;       % Tolerancia deseada
maxIter = 10000;     % maximo numero de iteraciones

size(A)


solexac=A\b;
[xSol, numIter, resvec] = cocg(A, b, x0, tol, maxIter);

fprintf('Solucion encontrada en %d iteraciones.\n', numIter);
fprintf('‖Xsol - solexac‖₂ = %.2e\n',norm(xSol - solexac));
