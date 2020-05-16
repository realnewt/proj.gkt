function y1 = idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x1)

x2=@(x1)(1-x1);

%gamma2= 1; %@(x1)(exp(-log(x1+W12.*x2(x1))+(x2(x1).*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));
%gamma1= 1;  %@(x1)(exp(-log(x2(x1)+W21.*x1)-(x1.*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));

P1sat=@(T)(exp(A1-B1./(T+C1))); P2sat=@(T)(exp(A2-B2/(T+C2)));

y1=@(T)(P1sat(T).*x1./P); y2=@(T)(P2sat(T).*x2(x1)./P);

sumy=@(T)(y1(T)+y2(T)-1);

T=fsolve(sumy,Tstart,optimset('Diagnostics','off', 'Display','off'));
y1=y1(T);
y2=y2(T);
end