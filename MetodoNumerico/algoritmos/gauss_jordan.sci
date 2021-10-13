function x = gauss_jordan(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Gauss-Jordan con pivoteo parcial.
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
//Normalizacion de E_k
for i =1:n
    a(i,:) = a(i,:)/a(i,i);
end
//Eliminacion de x_k en las ecuaciones por encima de la ecuacion k
for k=n-1:-1:1
    for i=n:-1:k+1
        a(k,:)= a(k,:)-a(k,i)*a(i,:);
    end
end
x = a(:,n+1);
endfunction
A = [0 2 3; 2 0 3; 8 16 -1]
b = [7 13 -3]'
disp("La solucion del sistema es: ")
disp(gauss_jordan(A,b));
