function res = findTbForDest(T, x1, W_Be_Ba, A1, A2, P)

x2 = 1 - x1;

gamma1 = exp(-log(x1)+W_Be_Ba(1)*x2)+x2.*((W_Be_Ba(1)./(x1+W_Be_Ba(1)*x2))-(W_Be_Ba(2)./(W_Be_Ba(2)*x1+x2)));
gamma2 = exp(-log(x2+W_Be_Ba(2)*x1)-x1.*((W_Be_Ba(1)./(x1+W_Be_Ba(1)*x2))-(W_Be_Ba(2)./(W_Be_Ba(2)*x1+x2))));

P01 = exp(A1(1) - A1(2)./(T + A1(3)));
P02 = exp(A2(1) - A2(2)./(T + A2(3)));

P1 = gamma1.*P01.*x1;
P2 = gamma2.*P02.*x2;

y1 = P1./P;
y2 = P2./P;

res = y1 + y2 - 1;
end

