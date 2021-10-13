function y = convergencia_iterativos(A,N)
    /*
    Toma una matriz A y su matriz del metodo iterativo a usar.
    Calcula si el metodo se pude asegurar convergente o no mediante si la matriz A
    es diagonalmente dominante, si alguna de las normas de I-N^{-1}*A o si su radio espectral es menor a 1.
    Devuelve 1 si se pude asegurar y 0 sino.
    */
    [n,m] = size(A);
    diagonal = 1;
    for i = 1:n;
        c =0;
        for j = 1:n
            if i ~= j then
                c = c+abs(A(i,j));
            end
            if c >= A(i,i) && diagonal then
                diagonal = 0;
                disp("La matriz no es diagonal dominante.");
                break;
            end
        end
        if diagonal == 0 then
            break;
        end
    end
    P = eye(n,n)-inv(N)*A;
    norm1 = norm(P,1);
    norm2 = norm(P,2);
    normi = norm(P,%inf);
    normmi = norm(P,-%inf);
    if norm1 < 1 || norm2 < 1 || normi < 1 && normmi < 1 then
        norma = 1;
    else
        norma = 0;
        disp("Las normas de P son mayores a 1");
    end
    e = max(abs(spec(P)));
    if e>1 then
        disp("El radio espectral es mayor que 1");
    end
    if diagonal || norma || e<1 then
        y = 1;
    else
            y = 0;
    end
endfunction
A = [1 -1 -1;0 2 4;1 -1 2];
N = [1,0,0;0,2,0;0,0,2];
if convergencia_iterativos(A,N) then
    disp("El sistema es convergente");
else
    disp("El metodo no se pude asegurar si es convergente");
end
A1 = [1,-1,0;-1,2,-1;0,-1,1.1];
N1 = [1,0,0;0,2,0;0,0,1.1];
if convergencia_iterativos(A1,N1) then
    disp("El metodo es convergente");
else
    disp("El metodo no se pude asegurar si es convergente");
end
