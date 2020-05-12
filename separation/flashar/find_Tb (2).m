function [T,y1]=find_Tb(P,T0,A1,A2,x1)
%function which finds T at which y1+y2=1 
x2=@(x1)(1-x1);

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3))));

y1=@(T)(P1sat(T).*x1./P); y2=@(T)(P2sat(T).*x2(x1)./P);

sumy=@(T)(y1(T)+y2(T)-1);

T=fsolve(sumy,T0,optimset('Diagnostics','off', 'Display','off'));
y1=y1(T);
y2=y2(T);
end