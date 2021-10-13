clc
clear
xdel(winsid())
// Entrada : x = numero real
// Salida : y = (e^x-1)/x - 1
function y = F(x)
    y = cos(x).*cosh(x)+1;
endfunction

//E1
x = linspace(0,8,5001);
y = F(x);
//plot(x,y,'b');
//legend("$F(x)$");
//Si vemos la grafica podes ver que hay ceros en 1.88 , 4.6 y 7.855

//E2
function y = pm(a,b)
    y = (a+b)/2;
endfunction
//f(x) = sin(x)-(x^2)/2
/*
a)
Veamos primero el rango entre [-0.5, 0.5].
f(-0.5)f(0.5) < 0
c = 0
b - c > epsilon => f(c) = 0, entonces es una raiz
Ahora veamos el rango [0.5, 3/2]
Ambos extremos tienen distinto signo
c = 1; b-c > epsilon => f(c) > 0 => a = c0
c1 = 1.25; b-c> E => f(c1)>0 => a = c1
c2 = 1.375; f(c2) > 0 => a=c2
c3 =1.4375; f(c3) < 0 => b=c3
c4 = 1.40625; f(c4) < 0 => b = c4
c5 = 1.390625; f(c5) >0 => a = c5
c6 = 1.3984375; b -c6 < epsilon => c6 es raiz
*/
function y = e2b(x)
    y = exp(-x)-x^4;
endfunction
/*
Veamos el rango [-1.5, -0.5]
c0= -1; f(c0)>0; b = c0
c1=-1.25; f(c1).0; b = c1
c2=-1.375; f(c2)>0; b = c2
c3=-1.4375; f(c3)< 0;a = c3
c4=-1.40625; f(c4)>0; b =c4
c5 = -1.421875; f(c5)>0; b=c5
c6=-1.4296875; b-c<epsilon => c6 raiz
veamos el rango [0.5, 1]
c0=0.75; f(c0)>0; a = 0.75
c1=0.875;f(c1)<0; b=0.875
c3=0.8125;f(c2)>0; a=0.8125
c4=0.84375;f(c4)<0; b=c4
c5=0.828125;f(c5)<0; b=c5
c6=0.8203125; b-c <10^-2; Raiz
*/
//E3
function y=e3(x)
    y =((x^2)/4)-sin(x);
endfunction
function y=ite3(x,xn)
    y = xn-e3(xn)*((xn-x)/(e3(xn)-e3(x)))
endfunction
/*
x0 = 1 y x1 = 2
x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0))) = 1.8670389
x3=1.9313546
x4=1.9338445
x5=1.9337536
x6=1.9337538
x7=1.9337538 por lo tanto converge
*/
//E4
/*
Sea x0 \in R un punto arbitrario. Sea la sucesion de aplicar reiteramente
a x0 la funcion coseno
x_{n+1}=cos(x_n), n \in N_0
Supongamos x_n \to a \in R. Luego, como cos(x) es continua, se tiene
a = cos(a)
Entonces a es un punto fijo de cos
Sea g(x)=cos(x), y sea x0 \in R arbitrario. Definimos
x_{n+1}= g(x_n), n \in N_0
Sea un intervalo [b,c] \subset R donde g([b,c])\subset [b,c] y {x_n}\subset[b,c].
Tomemos el intervalo [-1,1] se tiene
sup_{x\in[-1,1]} |g'(x)| = sup_{x\in[-1,1]} |-sin(x)| = sup_{x\in[-1,1]} |sin(x)| = sin(1) < 1
Como x_1=cos(x_0) \in [-1,1], por el teorema 2 aplicado a la sucesion {x_n:n>0}
tenemos que existe un unico punto a tal que x_n \to a cuando n \to infinito. Esto es a es un punto fijo
*/
//E5
function y=ge5_n(x0, n)
    y=x0
    for i=1:n
        y = 2^(y-1)
    end
endfunction
/*
Tomemos el intervalo [-inf, 2] se tiene que
-inf < x \leq 2 => -inf < 0 < g(x) <=2
Por lo tanto existe un x en [-inf,2] que vale. Si tomamos cualquier x0
menor o igual a 1, la iteracion converge y el limite es 1. Ya que si
calculamos el limite podemos ver q tiende a 1.
Si 1 < x0 < 2 entonces la iteracion converge en 1. Ya que si calculamos el limite seria como calcular el limte al infinito de 2^1/x.
Si x0=2 la iteracion converge a 2 ya que esta es solucion.
Si x0>2 entonces no converge ya que 2^k-1 > k
*/

//E6
/*
Se tiene que resolver x^2-5=0 es equivalente a buscar soluciones de x+c(x^2-5)=x para c>0.
Como g y g' son continuas, necesitamos que |g'(z)|<1 para asegurar la convergencia de x_n si x_0 es cercano a z.
|1-2 c sqrt(5)|<1 <=>-1<1-2 c sqrt(5)<1 <=>-2<-2 c sqrt(5)<0
1/sqrt(5)>c>0
Por lo tanto tomando 0<c<1/sqrt(5) podemos asegurar la convergencia de la sucesion para un valor adecuado de x0
*/

//E7
function y = ge7(d)
    y = (((2*%pi)/5)^2)/(9.8*tanh(4*d));
endfunction
function y = ge7_n(x0,n)
    y = x0;
    for i=1:n
        y = ge7(y)
    end
endfunction
// d = 0.2 con un solo digito de presicion => l = 10pi

function salida = newton(fun,x0,tol,iter)
    deff("y=f(x)","y="+fun);
    i = 0;
    x1 = x0 - f(x0)/numderivative(f,x0)
    while abs(x1-x0)>tol && i < iter
        i = i+1;
        x0 = x1;
        x1 = x0 - f(x0)/numderivative(f,x0)
    end
    if (abs(x1-x0)> tol) then disp('Se alcanzo el máximo de iteraciones'); end
    disp(i);
    salida = x1;
endfunction
s_e7 = "(((2*%pi)/5)^2)-(9.8*x*tanh(4*x))"
//d = 0.2250, entonces l = 27.9377

//E8
/*
### A ###
Tomemos el intervalo [0,1], vemos que g([0,1])\subset [0,1]. Prueba
Sea x \in [0,1], como g(x) es una funcion creciente,
g(0)\leq g(x) \leq g(1) => 1/3 \leq g(x) \leq e/3 => 0 \leq g(x) \leq 1
Luego g(x) \in [0,1]. De donde tenemos que g([0,1]) \subset [0,1].
sup_{x\in[0,1]} |g'(x)| =sup_{x\in[0,1]} e^x/3 = e/3 <1
Entonces, g verifica las hipotesis del teorema 2 en el intervalo y la sucesion
x_{n+1}=g(x_n)
converge a a_1. Por lo tanto la funcion es util para obtener una soluion.
### B ###
Consideremos g2(x) = (e^x-x)/2
e^x=3x <=> e^x-x = 2x <=> (e^x-x)/2 = x <=> g2(x)=x
Por lo tanto, una solucion sera un punto fijo de g2.
Sea x \in [0,1], como g2(x) es una funcion creciente para x \geq 0,
g2(0)\leq g2(x) \leq g2(1) => 1/2 \leq g(x) \leq (e-1)/2 => 0 \leq g(x) \leq 1
Luego g2(x) \in [0,1] => g2([0,1]) \subset [0,1]
g2'(x)= 1/2*(e^x-1)
sup_{x\in[0,1]} |g2'(x)|=sup_{x\in[0,1]} 1/2*(e^x-1) = 1/2*(e-1)<1
Entonces, g2 verifica las hipotesis del teorema 2 en el intervalo y la sucesion converge a a2. La funcion es util.
### c ###
Consideremos g3(x)= ln(3x)
e^x=3x <=> x = ln(3x) <=> x = g3(x)
Por lo tanto, una solucion sera un punto fijo de g3
Sea x \in [1.25,2], como g3 es creciente,
g3(1.25) \leq g3(x) \leq g3(2) => ln(3.75) \leq g3(x) \leq ln(6) => 1 \leq g3(x) \leq 2
Luego g3(x) \in [1.25,2] => g3([1.25,2]) \subset [1.25,2]
g3'(x)= 1/x
sup_{x\in[1.25,2]} |g3'(x)|=sup_{x\in[1.25,2]} 1/x =  0.8 < 1
Entonces, g3 verifica las hipotesis del teorema 2 en el intervalo y la sucesion converge a a3. La funcion es util.
### d ###
COnsideremos g4(x) = e^x-2x
e^x=3x <=> e^x-2x = x <=> g4(x) = x
Por lo tanto, una solucion sera un punto fijo de g4
Sea x \in [0.5.1], podemos ver que g4 es decreciente de [0.5,ln(2)]
g4(0.5) \geq g4(x) \geq g4(ln(2)) => 1 \geq g4(x) \geq 0.5
Podemos ver que g4 es creciente en [ln(2), 1]
g4(ln(2)) \geq g4(x) \geq g4(1) => 0 \geq g4(x) \geq e-2 \leq 1
Luego g4(x) \in [0.5,1] => g4([0.5,1])\subset [0.5,1]
g4'(x) = e^x-2
sup_{x\in[0.5,1]} |g4'(x)|=sup_{x\in[0.5,1]} |e^x-2|= e-2 < 1
Entonces, g4 verifica las hipotesis del teorema 2 en el intervalo y la sucesion converge a a3. La funcion es util.
*/
//E9
function y=fe9(x)
    f1 = 1+x(1)^2-x(2)^2+exp(x(1))*cos(x(2));
    f2 = 2*x(1)*x(2)+exp(x(1))*sin(x(2));
    y = [f1;f2];
endfunction
//deff("y=fe9(x)","y=[1+x(1)^2-x(2)^2+exp(x(1))*cos(x(2));2*x(1)*x(2)-4-x(2)^3]");
x0=[-1;4];
x1 = x0 - inv(numderivative(fe9,x0))*fe9(x0);
for i=1:4
    x0 = x1;
    x1 = x0 - inv(numderivative(fe9,x0))*fe9(x0);
end

//E10
function y=newtonmulti(f, x0, n)
    x1 = x0 - inv(numderivative(f,x0))*f(x0);
    for i = 1:(n-1)
        x0 = x1;
        x1 = x0 - inv(numderivative(f,x0))*f(x0);
    end
    y = x1;
endfunction
deff("y=fe10(x)","y=[x(1)^2+x(1)*x(2)^3-9;3*x(1)^2*x(2)-4-x(2)^3]")
xe10 = newtonmulti(fe10, [1.2;2.5], 10);
disp(string(fe10(xe10)));
xe10 = newtonmulti(fe10, [-2;2.5], 10);
disp(string(fe10(xe10)));
xe10 = newtonmulti(fe10, [-1.2;-2.5], 10);
disp(string(fe10(xe10)));
xe10 = newtonmulti(fe10, [2;-2.5], 10);
disp(string(fe10(xe10)));
//E11
deff("y=fe11(x)", "y=2*x(1)+3*x(2)^2+exp(2*x(1)^2+x(2)^2)");
deff("y=Dfe11(x)", "y=[2+4*x(1)*exp(2*x(1)^2+x(2)^2);6*x(2)+2*x(2)*exp(2*x(1)^2+x(2)^2)]"); //Df= grad f(x,y)
deff("y=Hfe11(x)","y=[16*(x(1)^2)*exp(2*x(1)^2+x(2)^2)+4*exp(2*x(1)^2+x(2)^2),8*x(1)*x(2)*exp(2*x(1)^2+x(2)^2);8*x(1)*x(2)*exp(2*x(1)^2+x(2)^2),4*x(2)^2*exp(2*x(1)^2+x(2)^2)+2*exp(2*x(1)^2+x(2)^2)+6]"); //Hf(x) = hessian de f(x,y)
x0 = [1;1];
x1 = x0 - inv(numderivative(Dfe11,x0))*Dfe11(x0);
while norm(x1-x0,2) > 10^-12
     x0 = x1;
     x1 = x0 - inv(numderivative(Dfe11,x0))*Dfe11(x0);
end
disp("El valor aproximado que satisface (i) es ");
disp(x1);
disp("El valor de Df(z) es ");
disp(Dfe11(x1));
disp("El valor de Hf(z) es ");
disp(numderivative(Dfe11, x1));
disp("Como sus autovalores son positivos, es una matriz definida positiva");
//E12
function y=newtonmulti_tol(f, x0, tol)
    x1 = x0 - inv(numderivative(f,x0))*f(x0);
    while norm(x1-x0,2) > tol
        x0 = x1;
        x1 = x0 - inv(numderivative(f,x0))*f(x0);
    end
    y = x1;
endfunction
/*
k(1)*exp(k(2))+k(3)-10 = 0
k(1)*exp(k(2)*2)+k(3)*2-12 = 0
k(1)*exp(k(2)*3)+k(3)*3-15 = 0
*/
deff("y=fe12(x)", "y=[x(1)*exp(x(2))+x(3)-10;x(1)*exp(x(2)*2)+x(3)*2-12;x(1)*exp(x(2)*3)+x(3)*3-15]");
k = newtonmulti_tol(fe12, [5;0.1;1], 10^(-2));
disp("El valor aproximado de k ");
disp(k)
disp("El valor de f(k) ");
disp(fe12(k));
//### B ###
deff("y=ge12(x)", "y="+string(k(1))+"*exp("+string(k(2))+"*x)+x*"+string(k(3))+"-500");
s_ge12 = string(k(1))+"*exp("+string(k(2))+"*x)+x*"+string(k(3))+"-500";
x = newton(s_ge12, 1, 10^(-2), 100);
disp("El valor del radio es ");
disp(x);
disp("El valor de f(d) es");
disp(ge12(x));
