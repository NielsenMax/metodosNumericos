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

function x = resol_cholesky(A,b)
    /*
    Utiliza el metodo de cholesky para calcular la solucion del
    sistema con la matriz de coeficientes A y el vector b.
    */
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
A = [16,-12,8;-12,18,-6;8,-6,8];
b = [76;-66;46];
disp("La solucion del sistema es: ")
disp(resol_cholesky(A,b));
