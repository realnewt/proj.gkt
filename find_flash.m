function res = find_flash(U, x1f, alfa, W12, W21,A1,B1,C1,A2,B2,C2,P)

T=U(1); x1 = U(2);
x2 = 1-x1;

gamma1 = exp(-log(x1+W12*x2)+x2.*((W12./(x1+W12*x2))-(W21./(W21*x1+x2))));
gamma2 = exp(-log(x2+W21*x1)-x1.*((W12./(x1+W12*x2))-(W21./(W21*x1+x2))));

P01 = 10.^(A1-B1./(T+C1));
P02 = 10.^(A2-B2./(T+C2));

P1 = gamma1.*P01.*x1;
P2 = gamma2.*P02.*x2;

y1 = P1./P;
y2 = P2./P;

res = [x1 + x1f - alfa
    y1+y2-1];

end