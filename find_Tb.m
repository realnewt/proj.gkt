function res = find_Tb(T,x1,A1,A2,P)

%Use T,x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P
%to calculate y1 and y2

x2 =1-x1;

P01 = 10.^(A1(1)-A1(2)./(T+A1(3)));
P02 = 10.^(A2(1)-A2(2)./(T+A2(3)));

P1 = P01.*x1;
P2 = P02.*x2;

y1 = P1./P;
y2 = P2./P;

res = y1+y2-1;
end

