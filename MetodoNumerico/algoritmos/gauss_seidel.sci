function x=resol_ti(L, b)
    /*
    Toma una matriz tringular inferior  nxn A y un vector n b
    Calcula x tal que A*x=b
    */
    [n,m] = size(L);
    b(1) = b(1)/L(1,1);
    for i=2:n
        b(i)= (b(i)-(L(i,1:(i-1))*b(1:(i-1))))/L(i,i);
    end
    x = b;
endfunction
//A2=[1,0,0;2,3,0;4,5,6];
//b2 = [1;8;32];
//disp("La solucion al sistema es ");
//disp(resol_ti(A2,b2));
function x = gauss_seidel(A,b,x0,tol,iter)
    /*
    Toma una matriz de coeficientes A nxn, dos vectores, un real y un natural.
    Calcula la solucion del sistema de ecuaciones determinado por la matriz de coeficientes a y el vector b.
    Partiendo de la solucion inicial x0 y aplicando el metodo de Gauss-seidel.
    La solucion estara dentro de la tolerancia tol y se calculara en menos de iter itercaiones.
    */
    [n,m] = size(A);
    [nb,mb] = size(b);
    if n<>m then
    error('gauss_seidel - La matriz A debe ser cuadrada');
    abort;
elseif m<>nb then
    error('gauss_seidel - dimensiones incompatibles entre A y b');
    abort;
end;
    N = tril(A);
    if det(N)==0 then
        error("gauss_seidel - La matriz era singular");
    end
    x1 = resol_ti(N,((N-A)*x0+b));
    i = 0;
    while norm(x1-x0)>tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = resol_ti(N,((N-A)*x0+b));
    end
    if (norm(x1-x0)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    x = x1;
endfunction
A = [10,1,2,3,4;1,9,-1,2,-3;2,-1,7,3,-5;3,2,3,12,-1;4,-3,-5,-1,15];
b = [12;-27;14;-17;12];
x=gauss_seidel(A,b,[1;1;1;1;1],10^-6,100)
disp("La solucion al sistema es: ");
disp(x);
