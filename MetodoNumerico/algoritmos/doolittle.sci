function [L,U]=doolittle(A)
    /*
    Utiliza el metodo de doolittle para calcular la factorizacion LU de la matriz A
    */
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
A= [16,-12,8,-16;-12,18,-6,9;8,-6,5,-10;-16,9,-10,46];
[L,U] = doolittle(A);
disp("La matriz triangular inferior es: ");
disp(L);
disp("La matriz triangular superior es: ");
disp(U);
function x =resol_doolittle(A, b)
        /*
    Utiliza el metodo de doolittle para calcular la solucion del
    sistema con la matriz de coeficientes A y el vector b.
    */
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
A= [16,-12,8,-16;-12,18,-6,9;8,-6,5,-10;-16,9,-10,46];
b = [-48;42;-29;156];
x = resol_doolittle(A,b);
disp("La solucion del sistema es: ");
disp(x);
