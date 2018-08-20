function g = Reinsch_0(Q,R,IwInv,Y,alpha)

T = R+alpha*(Q')*IwInv*Q;
L = chol(T,'lower');

gamma = (L')\(L\((Q')*Y));

g = Y-alpha*IwInv*Q*gamma;