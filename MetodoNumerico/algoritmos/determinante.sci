function d = determinante(A)
    /*
    Obtiene el determinante de la matriz A apartir de multiplicar
    los elementos de la diagonal principal de la matriz triangular
    superior resultante de la etapa de eliminacion progresiva
    del procedimiento de eliminacion gaussiana sin pivoteo.
    */
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
A = [1,2,3;3,-2,1;4,2,-1];
disp("El determinante de la matriz es: ")
disp(determinante(A));
