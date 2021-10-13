function [L,U,P]=factLUPP(A)
    /*
    Cacula la factorizacion LU de la matriz A, utilizando pivoteo parcial
    */
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
A = [2,1,1,0;4,3,3,1;8,7,9,5;6,7,9,8];
[L,U,P] = factLUPP(A);
disp("La matriz de permutacion es: ");
disp(P);
disp("La matriz triangular inferior es: ");
disp(L);
disp("La matriz triangular superior es: ");
disp(U);
