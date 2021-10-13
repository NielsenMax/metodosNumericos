function salida = newton(fun,x0,tol,iter)
    /*
    Toma un funcion en forma de string, dos reales y un natural.
    Encuentra una raiz para la funcion fun, aplicando el metodo de Newton partiendo de x0. 
    El resultado tendra una tolerancia de tol y el algoritmo hara como maximo iter iteraciones.
    */
    deff("y=f(x)","y="+fun);
    i = 0;
    x1 = x0 - f(x0)/numderivative(f,x0)
    while abs(x1-x0)>tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = x0 - f(x0)/numderivative(f,x0)
    end
    if (abs(x1-x0)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    salida = x1;
endfunction
fun = "(((2*%pi)/5)^2)-(9.8*x*tanh(4*x))";
deff("y=f(x)","y="+fun);
x0 = 0.2;
tol = 10^-6;
iter = 1000;
x = newton(fun,x0,tol,iter);
disp("El punto es x="+string(x)+". EL valor de es f(x)="+string(f(x)));


function y=newtonmulti(fun, x0, tol, iter)
    /*
    Toma un vector de n funciones en forma de string, un vector de n componentes reales,un real y un natural.
    Encuentra una raiz para el sistema de funciones fun, aplicando el metodo de Newton multivariable partiendo de x0. 
    El resultado tendra una tolerancia de tol y el algoritmo hara como maximo iter iteraciones.
    */
    deff("y=f(x)","y="+fun);
    i = 0;
    x1 = x0 - inv(numderivative(f,x0))*f(x0);
    while norm(x1-x0,2) > tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = x0 - inv(numderivative(f,x0))*f(x0);
    end
    if (norm(x1-x0,2)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    y = x1;
endfunction
fun ="[x(1)^2+x(1)*x(2)^3-9;3*x(1)^2*x(2)-4-x(2)^3]";
deff("y = g(x)","y="+fun);
x0 = [1.2; 2.5];
tol = 10^-6;
iter = 100;
x = newtonmulti(fun, x0, tol, iter);
disp("El punto es x=");
disp(x);
disp("EL valor es de f(x)=");
disp(g(x));
