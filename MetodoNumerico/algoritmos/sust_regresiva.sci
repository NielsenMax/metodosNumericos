
function x=resol_ts(A, b)
    /*
    Toma una matriz tringular superior  nxn A y un vector n b
    Calcula x tal que A*x=b
    */
    [n,m] = size(A);
    b(n) = b(n)/A(n,n);
     for i=(n-1):-1:1
            b(i)= (b(i)-(A(i,(i+1):n)*b((i+1):n)))/A(i,i);
    end
    x = b;
endfunction
A = [1,2,3;0,4,5;0,0,6];
b = [14;23;18];
disp("La solucion al sistema es ");
disp(resol_ts(A,b));

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
A2=[1,0,0;2,3,0;4,5,6];
b2 = [1;8;32];
disp("La solucion al sistema es ");
disp(resol_ti(A2,b2));
