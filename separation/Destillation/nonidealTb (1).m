function y1 = nonidealTb(P, Tstart, A1, A2, x1, W12)

x2 = @(x1) (1 - x1);

%Calculating non-ideality
gamma1=@(x1)(exp(-log(x1+W12(1).*x2(x1))+(x2(x1).*(W12(1)./(x1+W12(1).*x2(x1))-W12(2)./(W12(2).*x1+x2(x1))))));
gamma2=@(x1)(exp(-log(x2(x1)+W12(2).*x1)-(x1.*(W12(1)./(x1+W12(1).*x2(x1))-W12(2)./(W12(2).*x1+x2(x1))))));

%Saturation pressure
P1sat = @(T)(10.^(A1(1) - A1(2)./(T + A1(3))));
P2sat = @(T)(10.^(A2(1) - A2(2)./(T + A2(3))));

%Partial pressure
P1 = @(T) gamma1(x1).*(P1sat(T).*x1);
P2 = @(T) gamma2(x1).*(P2sat(T).*x2(x1));

%Mole fraction in the vapour phase
y1 = @(T) P1(T)./P;
y2 = @(T) P2(T)./P;

%Used to see if the sum of y1 and y2 is 1 as it should be
totY = @(T) (y1(T) + y2(T) - 1);

options = optimset('Display', 'off');

%Minimizing totY
T = fsolve(totY, Tstart, options);
y1 = y1(T);
end