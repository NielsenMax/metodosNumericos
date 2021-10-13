function x = jacobi(A,b,x0,tol,iter)
    /*
    Toma una matriz de coeficientes A nxn, dos vectores, un real y un natural.
    Calcula la solucion del sistema de ecuaciones determinado por la matriz de coeficientes a y el vector b.
    Partiendo de la solucion inicial x0 y aplicando el metodo de Jacobi.
    La solucion estara dentro de la tolerancia tol y se calculara en menos de iter itercaiones.
    */
    [n,m] = size(A);
    [nb,mb] = size(b);
    if n<>m then
    error('jacobi - La matriz A debe ser cuadrada');
    abort;
elseif m<>nb then
    error('jacobi - dimensiones incompatibles entre A y b');
    abort;
end;
    N = diag(diag(A));
    N_inv = eye(n,n)/N;
    if det(N)==0 then
        error("jacobi - La matriz era singular");
    end
    x1 = N_inv*((N-A)*x0+b);
    i = 0;
    while norm(x1-x0)>tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = N_inv*((N-A)*x0+b);
    end
    if (norm(x1-x0)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    x = x1;
endfunction
A = [10,1,2,3,4;1,9,-1,2,-3;2,-1,7,3,-5;3,2,3,12,-1;4,-3,-5,-1,15];
b = [12;-27;14;-17;12];
x=jacobi(A,b,zeros(5,1),10^-6,100)
disp("La solucion al sistema es: ");
disp(x);
