clc()
clear()
xdel()
//==================E1
//Funcion que toma un vector x y un real k, calcula la funcion usada en el polinomio de Lagrange usando los coeficientes del vector x y luego aplicandola en k.
function y=Lk(x,k)
    [Xn,Xm] = size(x);
    r = [x(1:k-1) x(k+1:Xm)];
    p = poly(r, "x", "roots");
    pk = horner(p, x(k));
    y = p/pk;
endfunction

//Funcion que toma dos vectores, uno a usar como nodos del polinomio y el otro
//como los valores de los nodos luego de aplicar la funcion a aproximar.
//La funcion calcula el polinomio de Lagrange correspondiente
function z = interpolLagrange(x,y)
    [Xn, Xm] = size(x)
    pol = 0;
    for k=1:Xm
        pol = pol + (Lk(x,k)*y(k));
    end
    z = pol
endfunction
/*
function p = newton(x, y, a)
    [n, m] = size(x);
    T(:,1) = x';
    T(:,2) = y';
    for j=3:(n+1)
        for i=(j-1):n
            T(i,j) = (T(i,(j-1))-T((i-1),(j-1)))/(T(i,1)-T((i-1),1));
        end
    end
    p=(a-T(n,1))*T(n,n+1);
    for i=(n-1):-1:2
        p = (a-T(i,1))*(T(i,i+1)+p);
    end
    p = T(1,2)+p;
endfunction
*/
//Funcion que aplica el metodo de diferencias divididas de Newton para aproximar una funcion
// Entrada: x,y = vectores puntos de interpolación (x,y)
// Salida: w = polinomio de diferencias divididas de Newton
function w = DD_Newton(x,y)
    w = 0
    s = poly(0,"x")
    n = length(x)
    for j=n:-1:2
        w = w*(s-x(j-1)) + DD(x(1:j),y(1:j))*(s-x(j-1))
    end
    w = w + y(1)
endfunction

//Funcion que calcula las diferencias divididas
// Entrada: x,y = vectores de puntos de interpolación (x,y)
// Salida: w = diferencias divididas en los vectores x e y
function w = DD(x,y)
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction
//lineal
x1 = [0, 0.2];
y1 = [1, 1.2214];
//cubica
x2  = [0, 0.2, 0.4, 0.6]
y2 = [1 1.2214 1.4918 1.8221];

lineal = interpolLagrange(x1,y1);
cubica = interpolLagrange(x2,y2);

rango = [-2:0.01:2];
//plot(rango, horner(lineal,rango), "r");
//plot(rango, horner(cubica, rango), "b");
//plot(rango, exp(rango), "g");
//a=gca();a.x_location = "origin";a.y_location = "origin";
//h1 = legend("lineal", "Cubico", "e^x");
/*
er_lineal(x) = (x-0)*(x-0.2)/2*exp(c_x)''
exp''(x) = exp(x)

er_lineal(1/3) = 0.0271423
*/
//=================4
// Funcion que calcula el error de interpolación en base a una cota para |f^(n)|
// Entrada: p = valor real, x = nodos de interpolación, cot = cota de |f^(n))|
// Salida: w = error de interpolación en x = p
function w = err(p,x,cot)
    n = length(x)
    w = cot/(factorial(n))
    for i=1:n do
        w = w*abs(p - x(i))
    end
endfunction
x4 = [2.0 2.1 2.2 2.3 2.4 2.5];
y4 = [0.2239 0.1666 0.1104 0.0555 0.0025 -0.0484];
p = DD_Newton(x4, y4);
w1 = horner(p, 2.15);
err1 = err(2.15, x4, 1);

w2 = horner(p, 2.35);
err2 = err(2.35, x4, 1);
//=============== E 7
// Funcion que calcula la aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo, tambien calcula su error
// Entrada: b = vectores 1xn
// Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
function [p,err] = MinCuad_pol(A,b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     err = A*w-b'
endfunction

//Funcion que calcula la matriz del método de mínimo cuadrados polinomial
// Entrada: x,y = vectores 1xn; n = grado de la aproximación
// Salida: A = matriz del método de mínimo cuadrados
function A = A_mc(x,n)
    m = length(x)
    A = ones(m,1)
    for j=2:(n+1) do
        A = [A,(x').^(j-1)]
    end
endfunction

// Método de Eliminación Gaussiana con pivoteo parcial
function [x,a] = gausselimPP(A,b)
[nA,mA] = size(A) 
[nb,mb] = size(b)
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
    temp = a(kpivot,:); a(kpivot,:) = a(k,:);
    a(k,:) = temp
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k)
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end
end
// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n)
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k)
    end
    x(i) = (a(i,n+1)-sumk)/a(i,i)
end
endfunction


x7= [0 0.15 0.31 0.5 0.6 0.75];
y7 = [1 1.004 1.31 1.117 1.223 1.422];
A = A_mc(x7,1);
deter = det((A')*A);
disp(deter);
[pe7_1, err7_1]=MinCuad_pol(A,y7);

A = A_mc(x7,2);
deter = det((A')*A);
disp(deter);
[pe7_2, err7_2]=MinCuad_pol(A,y7);

A = A_mc(x7,3);
deter = det((A')*A);
disp(deter);
[pe7_3, err7_3]=MinCuad_pol(A,y7);

// - Ejercicio 8 - //
x=[4,4.2,4.5,4.7,5.1,5.5,5.9,6.3,6.8,7.1]
y=[102.56,113.18,130.11,142.05,167.53,195.14,224.87,256.73,299.5,326.72]

disp("ítem a)")

disp("(#) n=1.")
A = A_mc(x,1)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 1 es:")
[p1,err1] = MinCuad_pol(A,y)
disp(p1)

disp("(#) n=2.")
A = A_mc(x,2)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 2 es:")
[p2,err2] = MinCuad_pol(A,y)
disp(p2)

disp("(#) n=3.")
A = A_mc(x,3)
deter = det((A')*A)
disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
disp("El polinomio de mínimos cuadrados de grado 3 es:")
[p3,err3] = MinCuad_pol(A,y)
disp(p3)


disp("(#) Analizamos los errores err = norm(Ax-y,2).")
disp("Para la aproximación lineal: "+string(norm(err1,2)))
disp("Para la aproximación cuadrática: "+string(norm(err2,2)))
disp("Para la aproximación cúbica: "+string(norm(err3,2)))
disp("Podemos decir que es mejor la aproximación cúbica en este caso.")

//B
rango = [3.5:0.01:7.5];
/*
plot(rango, horner(p1,rango), "r");
plot(rango, horner(p2, rango), "b");
plot(rango, horner(p3, rango), "g");
scatter(x,y);
*/
//============= E 9
//Funcion que calcula el polinomio de error de aproximacion
//Entrada: x =  vectores 1xn, cot = real que es la cota de |f^(n))|
//Salida: p = polinomio de error
function pe = poliErr(x,cot)
    n = length(x)
    pe = cot/(factorial(n));
    s = poly(0,"x");
    for i=1:n do
        pe = pe*abs(s - x(i));
    end
endfunction
function [p,e] = e9(a,b,n, cot)
    h = (b-a)/n;
    for i = 0:n
        x(i+1) = a + i * h;
    end
    y = 1./(1+x.^2);
    p = interpolLagrange(x,y);
    e = poliErr(x,cot);
endfunction
is = [2 4 6 10 14];
rango = [-5:0.01:5];
/*
|f''(x)|=|(6x^2-2)/(1+x^2)^3)|<= 2
*/
/*
[p9,er9]=e9(-5,5,is(1), 2);
plot(rango,horner(er9, rango),"r");
[p9,er9]=e9(-5,5,is(2), 2);
plot(rango,horner(er9, rango),"r");
[p9,er9]=e9(-5,5,is(3), 2);
plot(rango,horner(er9, rango),"r");
[p9,er9]=e9(-5,5,is(4), 2);
plot(rango,horner(er9, rango),"r");
[p9,er9]=e9(-5,5,is(5), 2);
plot(rango,horner(er9, rango),"r");
*/
/*
Podemos ver que al incrementar los nodos el error disminuye
*/
//========== E 10
//Funcion que calcula el polinomio n de Chebyshev
//Entrada: n = real que sera el grado del polinomio
//Salida: p = polinomio de Chebyshev correspondiente
function p = poliCheb(n)
    if n ==0 then
        p = 1;
    elseif n == 1
        p = poly(0,"x");
    else
        p = poly(0,"x");
        p = 2*p*poliCheb(n-1)-poliCheb(n-2);
    end
endfunction
x10 = roots(poliCheb(4));
y10 = exp(x10);
p10 = DD_Newton(x10,y10);
rango =[-1:0.01:1]
err10 = poliErr(x10, 2.72);
//plot(rango, horner(err10, rango), "r");
//=========== e11
x11 = roots(poliCheb(4));
x11 = 0.5*((%pi/2).*(x11+1));
y11 = cos(x10);
p11 = DD_Newton(x11,y11);
