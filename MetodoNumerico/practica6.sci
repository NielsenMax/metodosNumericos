clc()
clear()
xdel()
////==============E1
//Funcion que calcula los radios de los circulos de Gershgorin de una matriz
//Entrada: A = matriz nxn
//Salida: r = vector 1xn con los radios
function r = rGershgorin(A)
    [n,m] = size(A);
    for i=1:n
        r(i) = 0;
        for j=1:n
            if i ~= j
                r(i) = r(i)+abs(A(i,j));
            end
        end
    end
endfunction
A1_a = [1 0 0; -1 0 1;-1 -1 2];
//0 2 2
r1_a = rGershgorin(A1_a);
/*
Los autovalores se encuentran en {z: |z-1|<=0}U{z: |z-0|<=2}U{z: |z-2|<=2}
*/
A1_b = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2];
//0 0.2 0.2
r1_b = rGershgorin(A1_b);
disp(r1_b);
/*
Los autovalores se encuentran en {z: |z-1|<=0}U{z: |z-0|<=0.2}U{z: |z-2|<=0.2}
*/
A1_c = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2];
r1_c = rGershgorin(A1_c);
disp(r1_c);
/*
Los autovalores se encuentran en {z: |z-1|<=0}U{z: |z-0|<=0.5}U{z: |z-2|<=0.5}
*/
A1_d = [4 -1 0; -1 4 -1; -1 -1 4];
r1_d = rGershgorin(A1_d);
disp(r1_d);
/*
Los autovalores se encuentran en {z: |z-4|<=1}U{z: |z-4|<=2}U{z: |z-4|<=2}
*/
A1_e = [3 2 1; 2 3 0; 1 0 3];
r1_e = rGershgorin(A1_e);
disp(r1_e);
/*
Los autovalores se encuentran en {z: |z-3|<=3}U{z: |z-3|<=2}U{z: |z-3|<=1}
*/
A1_f = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75];
r1_f = rGershgorin(A1_f);
disp(r1_f);
/*
Los autovalores se encuentran en {z: |z-4.75|<=2.5}U{z: |z-4.75|<=3.5}U{z: |z-4.75|<=1.5}
*/
//---------B------
lamb1_a = spec(A1_a);
lamb1_b = spec(A1_b);
lamb1_c = spec(A1_c);
lamb1_d = spec(A1_d);
lamb1_e = spec(A1_e);
lamb1_f = spec(A1_f);
//========================E3
function A= genE3(n)
    A = [1 -1 0; -2 4 -2;0 -1 1];
    A(3,3) =  A(3,3)+0.1*n;
endfunction
A3_0=genE3(0);
A3_1=genE3(1);
A3_2=genE3(2);
A3_3=genE3(3);
A3_4=genE3(4);
A3_5=genE3(5);
A3_6=genE3(6);
A3_7=genE3(7);
A3_8=genE3(8);
A3_9=genE3(9);
A3_10=genE3(10);
p3_0 = poly(A3_0, 'x');
p3_1 = poly(A3_1, 'x');
p3_2 = poly(A3_2, 'x');
p3_3 = poly(A3_3, 'x');
p3_4 = poly(A3_4, 'x');
p3_5 = poly(A3_5, 'x');
p3_6 = poly(A3_6, 'x');
p3_7 = poly(A3_7, 'x');
p3_8 = poly(A3_8, 'x');
p3_9 = poly(A3_9, 'x');
p3_10 = poly(A3_10, 'x');
r3_0 = roots(p3_0);
r3_1 = roots(p3_1);
r3_2 = roots(p3_2);
r3_3 = roots(p3_3);
r3_4 = roots(p3_4);
r3_5 = roots(p3_5);
r3_6 = roots(p3_6);
r3_7 = roots(p3_7);
r3_8 = roots(p3_8);
r3_9 = roots(p3_9);
r3_10 = roots(p3_10);
lambda3_0 = spec(A3_0);
lambda3_1 = spec(A3_1);
lambda3_2 = spec(A3_2);
lambda3_3 = spec(A3_3);
lambda3_4 = spec(A3_4);
lambda3_5 = spec(A3_5);
lambda3_6 = spec(A3_6);
lambda3_7 = spec(A3_7);
lambda3_8 = spec(A3_8);
lambda3_9 = spec(A3_9);
lambda3_10 = spec(A3_10);
//======================E4
//Funcion que grafica un circulo centrado en x e y, con radio r
//Entrada: x,y = coordenadas del centro; r = radio del circulo;
function circ(x,y,r)
    xarc(x-r,y+r,r*2,r*2,0,360*64)
endfunction
//Funcion que grafica los circulos de Gershgorin de una matriz
//Entrada: A = matriz nxn
function Gers(A)
    [n,m] = size(A);
    centros = diag(A);
    radios = sum(abs(A),'c')-abs(centros);
    mx = round(min(centros-radios)-1);
    my = -round(max(abs(radios))+1);
    Mx = round(max(centros+radios)+1);
    My = -my;
    rectang = [mx,my,Mx,My];
    autov = spec(A);
    plot2d(real(autov), imag(autov), -1, "031", "", rectang);
    xgrid();
    for i =1:n
        circ(centros(i), 0,radios(i));
    end 
endfunction
/*
function Gers(A)
    [n,m] = size(A);
    for i=1:n
        r(i) = 0;
        for j=1:n
            if i ~= j
                r(i) = r(i)+abs(A(i,j));
            end
        end
    end
    maxx = max(diag(A));
    minx = min(diag(A));
    maxr = max(r);
    dim = max(abs(minx-maxr), abs(maxx+maxr));
    plot2d(1,1,strf="115",rect=[-dim,-dim,dim,dim]);
    a=gca(); // Handle on axes entity
    a.x_location = "origin"; 
    a.y_location = "origin"; 
    for i=1:n
        circ(A(i,i),0,r(i));
    end
endfunction
*/
//Funcion que grafica los circulos de Gershgorin junto con los autovalores de una matriz
//Entrada: A = matriz nxn
function CircGersValor(A)
    Gers(A);
    av = spec(A);
    x = real(av);
    y = imag(av);
    scatter(x,y);
endfunction
//=================E5
//Funcion que toma una matriz cuadrada, un vector, un real eps y un natural iter.
//Aplica el metodo de la potencia para calcular el mayor autovalor de la matriz.
//Realiza iter iteraciones como maximo y el algoritmo sigue hasta que el error sea menor a eps.
//Entrada: A = matriz nxn; z0 = vector nx1;
// eps = real, condicion de parada para el error;
// iter = natural, cota maxima de iteraciones
function [lambda,zn]=mPotencia(A,z0,eps,iter)
    autov = max(spec(A));
    lambda = 0;
    i = 1;
    w = A*z0;
    zn = w/norm(w,%inf);
    [m,j] = max(abs(w));
    lambda = m/z0(j);
    er1 = abs(autov-lambda);
    er2 = norm(z0-zn, %inf);
    er = max(er1, er2);
    while ((i<=iter) & (er>eps))
        z0 = zn;
        w = A*z0;
        zn = w/norm(w,%inf);
        [m,j] = max(abs(w));
        lambda = m/z0(j);
        er1 = abs(autov-lambda);
        er2 = norm(z0-zn, %inf);
        er = max(er1, er2);
        i = i+1;
    end
    disp("Se necesitaron ", i, " iteraciones");
endfunction
function [z,lambda]=metodoPotencia(A,z0,n,tol)
    w = A*z0;
    z1 = w/norm(w,'inf');
    i = 1;
    while i<n && norm(z1-z0) > tol
        i = i+1;
        z0 = z1;
        w = A*z0;
        z1 = w/norm(w,'inf');
    end
    for i =1:n
        if z0(i) ~= 0
            lambda = w(i)/z0(i);
            break;
        end
    end
    z = z1;
endfunction
A5_a = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6];
A5_b = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1];
[z_a,l_a]=metodoPotencia(A5_a, [1;0;0;0], 100, 10^-4);
[z_b,l_b]=metodoPotencia(A5_b, [1;0;0;0], 100, 10^-4);

