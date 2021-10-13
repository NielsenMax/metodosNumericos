clc
clear
xdel(winsid())
//Ejercicio 1
//Toma una matriz tringular superior  nxn A y un vector n b
//Calcula x tal que A*x=b
function x=resol_ts(A, b)
    [n,m] = size(A);
    b(n) = b(n)/A(n,n);
     for i=(n-1):-1:1
            b(i)= (b(i)-(A(i,(i+1):n)*b((i+1):n)))/A(i,i);
    end
    x = b;
endfunction
A = [1,2,3;0,4,5;0,0,6];
b = [14;23;18];
disp("La solucion al sistema es ");
disp(resol_ts(A,b));
//Toma una matriz tringular inferior  nxn A y un vector n b
//Calcula x tal que A*x=b
function x=resol_ti(A, b)
    [n,m] = size(A);
    b(1) = b(1)/A(1,1);
    for i=2:n
        for j= 1:(i-1)
            b(i)= b(i)-(A(i,j)*b(j));
        end
        b(i) = b(i)/A(i,i)
    end
    x = b;
endfunction
A2=[1,0,0;2,3,0;4,5,6];
b2 = [1;8;32];
disp("La solucion al sistema es ");
disp(resol_ti(A2,b2));

//Ejercicio 2
function [x,a] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

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

// Eliminación progresiva
n = nA;
for k=1:n-1
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
//### i ###
A2i = [1,1,0,3;2,1,-1,1;3,-1,-1,2;-1,2,3,-1];
b2i = [4;1;-3;4];
disp("Apartado i")
disp(gausselim(A2i,b2i));
//### ii ###
A2ii = [1,-1,2,-1;2,-2,3,-3;1,1,1,0;1,-1,4,3];
b2ii = [-8;-20;-2;4];
disp("Apartado ii")
[x,aa]=gausselim(A2ii,b2ii);
AA=[1,-1,2,-1;0,2,-1,1;0,0,-1,-1;0,0,0,-8];
disp(resol_ts(AA,b2ii));
//### iii ###
A2iii = [1,1,0,4;2,1,-1,1;4,-1,-2,2;3,-1,-1,2];
b2iii = [2;1;0;-3];
disp("Apartado iii")
disp(gausselim(A2iii,b2iii));
//### c ###
function [x,a] = gausselim_contadora(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  
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
c = 0;
// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            c = c + 1;
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
        c = c + 1;
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
disp("Se realizaron "+string(c)+" operaciones");
endfunction
//### d ###
function [x,a] = gausselim_sm(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  
[nA,mA] = size(A) 
[nb,mb] = size(b)
if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;
// Eliminación progresiva
n = nA;
for j = 1:(n-1)
    for i=(j+1):n
        mij = A(i,j)/A(j,j);
        A(i,j) = 0;
        A(i,(j+1):n)=A(i,(j+1):n)- mij*A(j,(j+1):n);
        b(i)= b(i)-mij*b(j);
    end
end
a = [A b];
b(n) = b(n)/A(n,n);
for i=(n-1):-1:1
    b(i)= (b(i)-(A(i,(i+1):n)*b((i+1):n)))/A(i,i);
end
x = b;
endfunction
//E3
function [x,a]=gausselim_multi(A,B)
    [nA,mA] = size(A) 
    [nB,mB] = size(B)
if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nB then
    error('gausselim - dimensiones incompatibles entre A y B');
    abort;
end;
// Eliminación progresiva
n = nA;
for j = 1:(n-1)
    for i=(j+1):n
        mij = A(i,j)/A(j,j);
        A(i,j) = 0;
        A(i,(j+1):n)=A(i,(j+1):n)- mij*A(j,(j+1):n);
        B(i,:)= B(i,:)-mij*B(j,:);
    end
end
a = [A B];
B(n,:) = B(n,:)./A(n,n);
for i=(n-1):-1:1
    B(i,:)= (B(i,:)-(A(i,(i+1):n)*B((i+1):n,:)))/A(i,i);
end
x = B;
endfunction
A = [1,2,3;3,-2,1;4,2,-1];
B = [14,9,-2;2,-5,2;5,19,12];
A_inv = gausselim_multi(A, eye(3,3));
//E4
function d = determinante(A)
        [n,m] = size(A)
if n<>m then
    error('determinante - La matriz A debe ser cuadrada');
    abort;
end;
d = 1;
for j = 1:(n-1)
    for i=(j+1):n
        mij = A(i,j)/A(j,j);
        A(i,j) = 0;
        A(i,(j+1):n)=A(i,(j+1):n)- mij*A(j,(j+1):n);
        d = d * A(j,j);
    end
end
d = d * A(n,n);
endfunction
A4 = [1,2,3;3,-2,1;4,2,-1];
//E5
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

// Ejemplo de aplicación
A2 = [0 2 3; 2 0 3; 8 16 -1]
b2 = [7 13 -3]'

[x2,a2] = gausselimPP(A2,b2)
disp(x2)
disp(a2)
//### i ###
disp("Apartado i")
disp(gausselimPP(A2i,b2i));
//### ii ###
disp("Apartado ii")
[x,aa]=gausselimPP(A2ii,b2ii);
disp(x);
//### iii ###
disp("Apartado iii")
disp(gausselim(A2iii,b2iii));
//-------------------E6---------------------
function [x,a]=gtd(A, b)
[nA,mA] = size(A) 
[nb,mb] = size(b)
if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;
n = nA;
c = 0;
for i = 2:n
    mij = A(i,(i-1))/A((i-1),(i-1));
    A(i,(i-1)) = 0;
    A(i,i:n)=A(i,i:n)- mij*A((i-1),i:n);
    b(i)= b(i)-mij*b(i-1);
    c = c+1;
end
a = A;
b(n) = b(n)/A(n,n);
c = c+1;
for i=(n-1):-1:1
    b(i)= (b(i)-(A(i,(i+1))*b(i+1)))/A(i,i);
    c = c+1;
end
x = b;
disp("Se necesitaron "+string(c)+" operaciones");
endfunction
Ae6=[1,2,0,0;3,4,5,0;0,6,7,8;0,0,9,10];
be6 = [5;26;65;67];
[xe6,ae6] = gtd(Ae6, be6);
disp("Solucion matriz tridiagonal ");
disp(ae6);
disp(xe6);
//----------------------------E7-------------------------
function [P,L,U]=factLUPP(A)
    U = A;
    [n,m] = size(A);
    if n<>m then
    error('factLUPP - La matriz A debe ser cuadrada');
    abort;
    end
    L = eye(n,n);
    P = eye(n,n);
    for k = 1:(n-1) do
        kpivot = k; umax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>umax then
                kpivot = i; umax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:n); U(kpivot,k:n) = U(k,k:n); U(k,k:n) = temp;
        temp = L(kpivot,1:(k-1)); L(kpivot,1:(k-1)) = L(k,1:(k-1)); L(k,1:(k-1)) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        for j=k+1:n do
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n)= U(j,k:n)-L(j,k)*U(k,k:m);
        end
    end
endfunction
Ae7 = [2,1,1,0;4,3,3,1;8,7,9,5;6,7,9,8];
[Pe7,Le7,Ue7] = factLUPP(Ae7);
disp(Pe7);
disp(Le7);
disp(Ue7);
//----------------------------E8----------------------
Ae8=[1.012,-2.132,3.104;-2.132,4.096,-7.013;3.104,-7.013,0.014];
Be8=[2.1756,4.0231,-2.1732,5.1967;-4.0231,6.0000,0,1.1973;-1.0000,5.2107,1.1111,0;6.0235,7.0000,0,4.1561];
[PA8,LA8,UA8] = factLUPP(Ae8);
[PB8,LB8,UB8] = factLUPP(Be8);
[LA8e,UA8e] = lu(Ae8);
[LB8e,UB8e] = lu(Be8);
// Son iguales solo que Le es la transpuesta

//---------------------------E9----------------------
A9 = [1,2,-2,1;4,5,-7,6;5,25,-15,-3;6,-12,-6,22];
b9= [1;2;0;1];
[LA9, UA9, PA9] = lu(A9);
x9 = resol_ts(UA9, resol_ti(LA9, PA9*b9));
//### B ####
b9b=[2;2;1;0];
x9b = resol_ts(UA9, resol_ti(LA9, PA9*b9b));

//---------------------------E10---------------------
function [L,U]=doolittle(A)
    [n,m] = size(A);
    if n<>m then
    error('doolittle - La matriz A debe ser cuadrada');
    abort;
    end
    U(1,:) = A(1,:);
    L(:,1)= A(:,1)/A(1,1);
    for i = 2:n
        U(i,:) = A(i,:)-L(i,1:(i-1))*U(1:(i-1),:);
        L(:,i) = (A(:, i)-L(:,1:(i-1))*U(1:(i-1),i))/U(i,i);
    end
endfunction
function x =resol_doolittle(A, b)
    [L,U] = doolittle(A);
    [n,m] = size(A);
    for i=2:n
        b(i)= b(i)-(L(i,1:(i-1))*b(1:(i-1)));
    end
    b(n) = b(n)/U(n,n);
    for i=(n-1):-1:1
            b(i)= (b(i)-(U(i,(i+1):n)*b((i+1):n)))/U(i,i);
    end
    x = b;
endfunction
//-------------------------E11----------------
function [U, ind] = Cholesky(A)
eps = 1.0e-8
n = size(A,1)
U = zeros(n,n)
for k = 1:n
    if k==1 then
            t = A(k,k)
    else 
            t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    end

    if t <= eps
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    U(k,k)= sqrt(t)
    for j = k+1:n
        if k==1 then 
                    U(k,j) = A(k,j)/U(k,k)
        else 
                    U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
        end
    end
end
ind = 1
endfunction
A11 = [16,-12,8,-16;-12,18,-6,9;8,-6,5,-10;-16,9,-10,46];
B11 = [4,1,1;8,2,2;1,2,3];
C11 = [1,2;2,4];
[A11C, indA11] = Cholesky(A11);
[B11C, indB11] = Cholesky(B11);
[C11C, indC11] = Cholesky(C11);
//------------------E12-------------
function x = resol_cholesky(A,b)
    [n,m] = size(A);
    if n<>m then
    error('resol_cholesky - La matriz A debe ser cuadrada');
    abort;
    end
    [U, ind] = Cholesky(A);
    if ind == 0 then
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    L = U'
    b(1) = b(1)/L(1,1);
    for i=2:n
        b(i)= (b(i)-(L(i,1:(i-1))*b(1:(i-1))))/L(i,i);
    end
    b(n) = b(n)/U(n,n);
    for i=(n-1):-1:1
            b(i)= (b(i)-(U(i,(i+1):n)*b((i+1):n)))/U(i,i);
    end
    x = b;
endfunction
A12 = [16,-12,8;-12,18,-6;8,-6,8];
b12 = [76;-66;46];
disp(resol_cholesky(A12,b12));
