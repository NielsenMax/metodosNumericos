function [Q,R]=fact_QR(A)
    /*
    Toma una matriz A de nxm y calcula su factorizacion QR.
    */
    [n,m] = size(A)
    if m > n  then
        error('fact_qr - La matriz A nxm debe ser n <= m');
        abort;
    end
    Q(:,1) = A(:,1)/norm(A(:,1));
    R(1,1) = norm(A(:,1));
    for k = 2:m
        sumk = zeros(n,1);
        for i = 1:k-1
            sumk = sumk + (A(:,k)'*Q(:,i))*Q(:,i)
        end
        R(k,k)=norm(A(:,k)-sumk);
        Q(:,k)=(A(:,k)-sumk)/norm(A(:,k)-sumk);
        for i=1:m
        for j=i+1:m
            R(i,j) = A(:,j)'*Q(:,i);
        end
    end
    end
endfunction
A = [1 0 5;1 3 1;2 -1 7];
[Q,R] = fact_QR(A);
disp("La matriz Q resulta: ");
disp(Q);
disp("La matriz R resulta: ");
disp(R);
disp(Q*R);
