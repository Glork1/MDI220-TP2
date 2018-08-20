% cd /home/roueff/svn/perso/trunk/stat1/tex/tp/non-param-reg/script/
clear all;
clc;
load Motorcycledata.txt
time = Motorcycledata(:,2);
Y = Motorcycledata(:,3);

close all;
plot(time,Y,'o');
the_legend{1} = 'données';
legend(the_legend);

% Reduction aux instants distincts

[alltime,ind,indrev] = unique(time);

for k = 1:length(alltime)
    Ym(k) = mean(Y(indrev==k));
    w(k) = length(find(indrev==k));
end

hold on
plot(alltime,Ym,'xr');
the_legend{2} = 'données après réduction des multiples';
legend(the_legend);

% Régresssion polynomiale

i = 1;
for d_p = [4 8 16]
    pol = polyfit(time,Y,d_p);
    pol_int = polyval(pol,time);
    hold on
    plot(time,pol_int,'.','Color',rand(1,3));
    the_legend{2+i} = ['polynôme de degré ' num2str(d_p)];
    i = i+1;
end
legend(the_legend);

% Régresssion spline

n = length(Ym);
p = length(alltime)-1; % doit etre égal à n-1!
h = diff(alltime);

Q = zeros(p+1,p-1);
sQ1 = zeros(3,p-1);
for j=1:(p-1)
    Q(j:(j+2),j) = [1/h(j);-1/h(j)-1/h(j+1);1/h(j+1)];
    sQ1(:,j) = [1/h(j);-1/h(j)-1/h(j+1);1/h(j+1)];
end

R = diag(h(1:end-1)+h(2:end))/3;
for j=1:p-2;
    R(j,j+1) = h(j+1);
    R(j+1,j) = h(j+1);
end

Iw = diag(w);
w_inv = 1./w;
IwInv = diag(w_inv);
sIwInv = sparse(1:n,1:n,w_inv,n,n);
sQ = sparse(Q);
sR = sparse(R);
K = sQ*(sR\eye(p-1))*(sQ');

for alpha = [10 20 30]
	g = (Iw+alpha*K)^(-1)*Iw*Ym(:);
    plot(alltime,g,'-.','Color',rand(1,3));
    the_legend{2+i} = ['spline \alpha=' num2str(alpha)];
    i = i+1;
end
legend(the_legend);    
    
% CV

alphas = 1:100;
n_alphas = length(alphas);
cv = zeros(1,n_alphas);

tic;

% [V,D] = eig(IwInv*K);                                         % 2
% V_inv = V\eye(p+1);                                           % 2
% d_D = diag(D);                                                  % 2
% f_A = @(beta) V*diag((1+beta*d_D).^(-1))*V_inv;                 % 2

for i = 1:n_alphas
    alpha = alphas(i);
%     A = ((Iw+alpha*K)^(-1))*Iw;                               % 1 pire
%     A = f_A(alpha);                                           % 2
%     g = A*Ym(:);                                              % 1 et 2
%     cv(i) = w*(((Ym(:)-g)./(1-diag(A))).^2);                  % 1 et 2
    [d g] = Reinsch(sQ,sQ1,sR,sIwInv,w_inv,Ym(:),alpha,1);    % 5 optimal
    cv(i) = w*(((Ym(:)-g)./d).^2);                            % 5 optimal
end

[aa,alpha_min_cv] = min(cv);
time_cv = toc;

close all;
plot(cv);

alpha = alpha_min_cv;
g_cv = (Iw+alpha*K)^(-1)*Iw*Ym(:);
plot(alltime,g_cv,'k',time,Y,'or');

% GCV

gcv = zeros(1,n_alphas);

tic;

d = eig(IwInv*K);                                               % 3, 4 et 5
den = @(beta) n*(1-mean((1+beta*d).^(-1)))^2;                   % 3, 4 et 5

% [V,D] = eig(IwInv*K);                                           % 2
% V_inv = V\eye(p+1);                                             % 2
% d = diag(D);                                                    % 2
% f_A = @(beta) V*diag((1+beta*d).^(-1))*V_inv;                   % 2

for i = 1:n_alphas
    alpha = alphas(i);
%     A = ((Iw+alpha*K)^(-1))*Iw;                               % 1 pire
%     A = f_A(alpha);                                           % 2
%     g = A*Ym(:);                                              % 1 et 2
%     g = Reinsch_0(Q,R,IwInv,Ym(:),alpha);                     % 3 Reinsch
%     g = Reinsch_0(sQ,sR,sIwInv,Ym(:),alpha);                     % 4 Reinsch sparse
%     den = (1-trace(A)/n)^2;                                   % 1 et 2
%     gcv(i) = w*((Ym(:)-g).^2)/den;                            % 1 et 2 
    [d_1 g] = Reinsch(sQ,sQ1,sR,sIwInv,w_inv,Ym(:),alpha,0);  % 5 optimal
    gcv(i) = w*((Ym(:)-g).^2)/den(alpha);                     % 3, 4 et 5
end

[aa,alpha_min_gcv] = min(gcv);
time_gcv = toc;

close all;
plot(gcv);

alpha = alpha_min_gcv;

g_gcv = (Iw+alpha*K)^(-1)*Iw*Ym(:);
plot(alltime,g_gcv,'k',time,Y,'or');

close all;
plot(time,Y,'o',alltime,g_cv,alltime,g_gcv);
legend('données',['CV, temps de calcul = ' num2str(time_cv) ' secondes'],['GCV, temps de calcul = ' num2str(time_gcv) ' secondes']);