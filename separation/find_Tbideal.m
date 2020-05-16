function res = find_Tbideal(T,x1,x2,A1,B1,C1,A2,B2,C2,P)

%Use T,x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P
%to calculate y1 and y2

P01 = exp(A1-B1./(T+C1));
P02 = exp(A2-B2./(T+C2));

P1 = P01.*x1; 
P2 = P02.*x2;

y1 = P1./P;
y2 = P2./P;

res = y1+y2-1;
end