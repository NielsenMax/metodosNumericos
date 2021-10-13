function w = optw_SOR(A)
    /*
    Dada una matriz definida positiva y tridiagonal, calcula el factor de escala optimo para el metodo SOR.
    */
            [n,m] = size(A);
if n<>m then
    error('gauss_seidel - La matriz A debe ser cuadrada');
    abort;
end;
D = diag(diag(A));
Tj = eye(n,n)-inv(D)*A;
eigen = spec(Tj);
maxe = abs(eigen(1));
for i=2:size(eigen,1)
    if abs(eigen(i)) > maxe then
        maxe= abs(eigen(i));
    end
end
w = 2/(1+sqrt(1-maxe^2));
endfunction
A = [4,3,0;3,4,-1;0,-1,4];
w = optw_SOR(A);
disp("El valor optimo de w es: ");
disp(w);


function [x]=gauss_seidel_SOR(A, b, x0, w, tol, iter)
    /*
    Dada una matriz de coeficientes A y el vector b, calcula la solucion
    del sistema Ax=b utilizando el metodo SOR, siendo w el factor de escala.
    La solucion tendra una tolerancia tol y se calculara usando menos de iter iteraciones.
    */
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
Tw = inv(D+w*L)*((1-w)*D-w*U);
cw = w*inv(D+w*L)*b;
x1 = Tw*x0+cw;
i = 0;
while norm(x1-x0) > tol && i < iter do
    i = i+1;
    x0 = x1;
    x1 = Tw*x0+cw;
end
if (norm(x1-x0)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
x = x1;
endfunction
A = [4,3,0;3,4,-1;0,-1,4];
b = [24;30;-24];
disp("La solucion del sistema es: ");
[x] = gauss_seidel_SOR(A,b,zeros(3,1),optw_SOR(A),10^-7,100);
disp(x);
