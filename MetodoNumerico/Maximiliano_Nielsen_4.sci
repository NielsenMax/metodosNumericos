clc
clear
xdel(winsid())
// Entrada : x = numero real
// Salida : y = (e^x-1)/x - 1
function y = F(x)
    y = (((%e.^x)-1)./x)-1;
endfunction

//a
x = linspace(-10e-8,10e-8,2001);
y = F(x);
plot(x,y,'b');
legend("$F(x)$");
//Imprimira una grafica con mucho ruido
//El resultado de la grafica es erroneo ya que es una propagacion de errores desde que el calculo e^x es
//una aproximacion. Tambien se añade el hecho de que e^x para numeros cercanos a cero es igual a 1
//lo que hace una resta de numeros similares que produce mayor error.
//Luego se vuelve a hacer otra resta de numeros similares cuando se realiza (e^x-1)/x - 1 ya que
//para x cercanos a 0 la primer parte de la resta tiende a uno.

//b
//Si analizamos la expresion de forma matematica la funcion tiende a 0 cuando x tiende a cero.
//Si analizamos el resultado de la grafica de scilab pareciera que tambien tendemos a 0 aunque es dificil de decidir.

//c
//Como se dijo en el apartado A los resultados cercanos a cero son erroneos dadas las restas de numeros similares que tenemos.
//Por ejemplo si miramos el error relativo de x = 2e-9, comparado con el valor "real" calculado por wolframalpha, podemos ver
//que el error relativo es de 29.281932.
disp("F(2e-9) ="+string(F(2e-9)))
disp("Error comparado con wolframalpha = "+string(abs(F(2e-9)-1e-9)/1e-9))
//Salida:
//F(2e-9) =-2.828D-08
//Error comparado con wolframalpha = 29.281932
