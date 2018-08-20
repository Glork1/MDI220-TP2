function [d g] = Reinsch(sQ,sQ1,sR,sIwInv,w_inv,Y,alpha,CV)

n = length(w_inv)-2;

sT = sR+alpha*(sQ')*sIwInv*sQ;
L = chol(sT,'lower');

gamma = (L')\(L\((sQ')*Y));

g = Y-alpha*sIwInv*sQ*gamma;

if (CV)
    L_diag_inv = 1./diag(L);
    d_aux = 1./diag(L);
    sTInv = zeros(n,3);
    sTInv(n,1) = L_diag_inv(n)*d_aux(n);
    sTInv(n-1,2) = -L_diag_inv(n-1)*L(n,n-1)*sTInv(n,1);   
    sTInv(n-1,1) = L_diag_inv(n-1)*(d_aux(n-1)-L(n,n-1)*sTInv(n-1,2));

    for i = (n-2):-1:1
        sTInv(i,3) = -L_diag_inv(i)*(L(i+1,i)*sTInv(i+1,2)+L(i+2,i)*sTInv(i+2,1));
        sTInv(i,2) = -L_diag_inv(i)*(L(i+1,i)*sTInv(i+1,1)+L(i+2,i)*sTInv(i+1,2));
        sTInv(i,1) = L_diag_inv(i)*(d_aux(i)-L(i+1,i)*sTInv(i,2)-L(i+2,i)*sTInv(i,3));
    end
    d = zeros(n+2,1);
    d(1) = w_inv(1)*((sQ1(1,1))^2)*sTInv(1,1);
    d(2) = w_inv(2)*(((sQ1(2,1))^2)*sTInv(1,1)+((sQ1(1,2))^2)*sTInv(2,1)+2*sQ1(2,1)*sQ1(1,2)*sTInv(1,2));
    d(n+1) = w_inv(n+1)*(((sQ1(3,n-1))^2)*sTInv(n-1,1)+((sQ1(2,n))^2)*sTInv(n,1)+2*sQ1(3,n-1)*sQ1(2,n)*sTInv(n-1,2));
    d(n+2) = w_inv(n+1)*((sQ1(3,n))^2)*sTInv(n,1);
    for i = 3:n
        d(i) = w_inv(i)*(((sQ1(3,i-2))^2)*sTInv(i-2,1)+((sQ1(2,i-1))^2)*sTInv(i-1,1)+((sQ1(1,i))^2)*sTInv(i,1)+2*sQ1(3,i-2)*sQ1(2,i-1)*sTInv(i-2,2)+2*sQ1(3,i-2)*sQ1(1,i)*sTInv(i-2,3)+2*sQ1(2,i-1)*sQ1(1,i)*sTInv(i-1,2));
    end
    d = alpha*d;
else
    d = 0;
end