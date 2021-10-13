function [x,a] = gausselim_sm(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo. Aprovechando los calculos de submatrices
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
A = [1,1,0,3;2,1,-1,1;3,-1,-1,2;-1,2,3,-1];
b = [4;1;-3;4];
disp("Solucion del sistema: ");
disp(gausselim_sm(A,b));

function [x,a]=gausselim_multi(A,B)
    /*
     Esta función obtiene la solución del sistema de ecuaciones lineales A*x=B, 
    dada la matriz de coeficientes A nxn y la matriz B nxm.
    La función implementa el método de Eliminación Gaussiana sin pivoteo.
    */
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
[x, a] = gausselim_multi(A,B);
disp("La matriz resultante de la eliminacion: ");
disp(a);
disp("La solucion del sistema es: ");
disp(x);
