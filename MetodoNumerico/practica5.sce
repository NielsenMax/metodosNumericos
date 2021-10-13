clc
clear
xdel(winsid())
//------------E1-----------------
A1 = [1 -1 -1;0 2 4;1 -1 2];
//A1 = [2,0,4;-1,1,-1;-1,1,2];
N1 = [1,0,0;0,2,0;0,0,2];
//N1 = [2,0,0;0,1,0;0,0,2];
c1 = norm(eye(3,3)-inv(N1)*A1,2);
e1 = max(abs(spec(eye(3,3)-inv(N1)*A1)));
//A no es diagonal dominante
//||I-N^[-1]*A|| > 1
//El radio esprectral de la matriz I-N^[-1]*A es mayor que uno.
//Por lo tanto no se puede garantizar la convergencia.

A1ii = [1,-1,0;-1,2,-1;0,-1,1.1];
N1ii = [1,0,0;0,2,0;0,0,1.1];
c1ii = norm(eye(3,3)-inv(N1ii)*A1ii,2);
e1ii = max(abs(spec(eye(3,3)-inv(N1ii)*A1ii)));
//A no es diagonal dominante
//||I-N^[-1]*A|| > 1
//El radio esprectral de la matriz I-N^[-1]*A es menor que uno.
//Por lo tanto se puede garantizar la convergencia.
//### B ###
N1b = [1,0,0;0,2,0;1,-1,2];
c1b = norm(eye(3,3)-inv(N1b)*A1);
e1b = max(abs(spec(eye(3,3)-inv(N1b)*A1)));
//A no es diagonal dominante
//||I-N^[-1]*A|| > 1
//El radio esprectral de la matriz I-N^[-1]*A es mayor que uno.
//Por lo tanto no se puede garantizar la convergencia.
N1bii= [1,0,0;-1,2,0;0,-1,1.1];
c1bii = norm(eye(3,3)-inv(N1bii)*A1ii);
e1bii = max(abs(spec(eye(3,3)-inv(N1bii)*A1ii)));
//A no es diagonal dominante
//||I-N^[-1]*A|| > 1
//El radio esprectral de la matriz I-N^[-1]*A es menor que uno.
//Por lo tanto se puede garantizar la convergencia.
//### c ###
function x=resol_ti(L, b)
    [n,m] = size(L);
    b(1) = b(1)/L(1,1);
    for i=2:n
        b(i)= (b(i)-(L(i,1:(i-1))*b(1:(i-1))))/L(i,i);
    end
    x = b;
endfunction
x0 = [0;0;0];
x1 = resol_ti(N1b,(N1b-A1)*x0+[0.375;0;0]);
i = 0
while norm(x1-x0) > 10^-2 && i<25 do
    i=i+1;
    x0 = x1;
    x1 = resol_ti(N1b,(N1b-A1)*x0+[0.375;0;0]);
end
disp(x1);
x0 = resol_ti(N1bii,(N1bii-A1ii)*[0;0;0]+[0;1;0]);
x1 = resol_ti(N1bii,(N1bii-A1ii)*x0+[0;1;0]);
while norm(x1-x0) > 10^-2 do
    x0 = x1;
    x1 = resol_ti(N1bii,(N1bii-A1ii)*x0+[0;1;0]);
end
disp(x1);
//------------------------E2------------------
A2 = [10,1,2,3,4;1,9,-1,2,-3;2,-1,7,3,-5;3,2,3,12,-1;4,-3,-5,-1,15];
N2j= [10,0,0,0,0;0,9,0,0,0;0,0,7,0,0;0,0,0,12,0;0,0,0,0,15];
b2 = [12;-27;14;-17;12];
c2 = norm(eye(5,5)-inv(N2j)*A2);
x0 = resol_ti(N2j,(N2j-A2)*[0;0;0;0;0]+b2);
x1 = resol_ti(N2j,(N2j-A2)*x0+b2);
cont =0
while norm(x1-x0) > 10^-6 do
    x0 = x1;
    x1 = resol_ti(N2j,(N2j-A2)*x0+b2);
    cont = cont +1;
end
disp(A2*x1);
disp(cont)
N2gs = [10,0,0,0,0;1,9,0,0,0;2,-1,7,0,0;3,2,3,12,0;4,-3,-5,-1,15];
x0 = resol_ti(N2j,(N2gs-A2)*[0;0;0;0;0]+b2);
x1 = resol_ti(N2j,(N2gs-A2)*x0+b2);
cont =0
while norm(x1-x0) > 10^-6 do
    x0 = x1;
    x1 = resol_ti(N2gs,(N2gs-A2)*x0+b2);
    cont = cont +1;
end
disp(A2*x1);
disp(cont)
//-----------------------E3--------------
n = 5;
A3 = 2*eye(n,n)- [[zeros(1,n-1);eye(n-1,n-1)] zeros(n,1)] - [zeros(n,1) [eye(n-1,n-1);zeros(1,n-1);]]
N3 = tril(A3);
P3 = eye(n,n)-inv(N3)*A3;
//----------------------E4---------------
function [x,a] = gausselimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k);
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
endfunction
function [x,c] = gauss_seidel(A,b,x0, tol)
    [n,m] = size(A);
    [nb,mb] = size(b)
if n<>m then
    error('gauss_seidel - La matriz A debe ser cuadrada');
    abort;
elseif n<>nb then
    error('gauss_seidel - dimensiones incompatibles entre A y b');
    abort;
end;
N=tril(A);
if det(N) == 0 then
    error('gauss_seidel - La matriz N debe resultar no singular');
    abort;
end
x1 = resol_ti(N,(N-A)*x0+b);
c =1;
while norm(x1-x0) > tol do
    x0 = x1;
    x1 = resol_ti(N,(N-A)*x0+b);
    c = c+1;
end
x = x1;
endfunction
function [A,b]=genE4(n)
    A = 8*eye(n,n) + 2*diag(ones(n-1,1),1)+2*diag(ones(n-1,1),-1)
+ diag(ones(n-3,1),3) + diag(ones(n-3,1),-3);
    b = ones(n,1);
endfunction
n = 100;
[A4,b4] = genE4(n);
tic();
x = gausselimPP(A4,b4);
t = toc();
disp("Eliminiacion gauss: ");
disp(t);
tic();
x = gauss_seidel(A4,b4,zeros(n,1),10^-6);
t = toc();
disp("Gauss Seidel 10^-6");
disp(t);
tic();
x = gauss_seidel(A4,b4,zeros(n,1),10^-11);
t = toc();
disp("Gauss Seidel 10^-11");
disp(t);
n = 500;
[A4,b4] = genE4(n);
tic();
//x = gausselimPP(A4,b4);
t = toc();
disp("500: Eliminiacion gauss: ");
disp(t);
tic();
x = gauss_seidel(A4,b4,zeros(n,1),10^-6);
t = toc();
disp("500: Gauss Seidel 10^-6");
disp(t);
tic();
x = gauss_seidel(A4,b4,zeros(n,1),10^-11);
t = toc();
disp("500: Gauss Seidel 10^-11");
disp(t);
//-------------------E5------------------
function [x,c]=gauss_seidel_SOR(A, b, x0, tol)
        [n,m] = size(A);
    [nb,mb] = size(b)
if n<>m then
    error('gauss_seidel - La matriz A debe ser cuadrada');
    abort;
elseif n<>nb then
    error('gauss_seidel - dimensiones incompatibles entre A y b');
    abort;
end;
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
Tj = eye(n,n)-inv(D)*A;
eigen = spec(Tj);
maxe = abs(eigen(1));
for i=2:size(eigen,1)
    if abs(eigen(i)) > maxe then
        maxe= abs(eigen(i));
    end
end
w = 2/(1+sqrt(1-maxe^2));
Tw = inv(D+w*L)*((1-w)*D-w*U);
cw = w*inv(D+w*L)*b;
x1 = Tw*x0+cw;
c = 1;
while norm(x1-x0) > tol do
    x0 = x1;
    x1 = Tw*x0+cw;
    c = c +1;
end
x = x1;
endfunction
//### a ###
A5 = [4,3,0;3,4,-1;0,-1,4];
b5 = [24;30;-24];
[x,c] = gauss_seidel(A5,b5,zeros(3,1),10^-7);
disp(x);
disp("Iteraciones necesarias: "+string(c));
//### b ###
[x,c] = gauss_seidel_SOR(A5,b5,zeros(3,1),10^-7);
disp(x);
disp("Iteraciones necesarias: "+string(c));
