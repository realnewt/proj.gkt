%%
%Binary dest. of buthene (1) and buthane (2)

clear, clc, clf

%Wilson parametrar
W12 = [0.48584, 1.64637];

%Antoine constants for log10, K, bar
A1 = [3.6470, 799.055, -46.615];
A2 = [4.3281, 1132.108, 0.918];


P = 1; %Bar

%Saturation pressure
Tb1 = A1(2)/(A1(1) - log10(P)) - A1(3);
Tb2 = A2(2)/(A2(1) - log10(P)) - A2(3);

x1 = linspace(0, 1);

options = optimset('Display', 'off');

%Calc of composition
for i = 1:length(x1)
    
    x2=@(x1)(1-x1(i));
    
    gamma1 = @(x1) (exp(-log(x1(i) + W12(1).*x2(x1)) + (x2(x1).*(W12(1)./(x1(i) + W12(1).*x2(x1)) - W12(2)./(W12(2).*x1(i) + x2(x1))))));
    gamma2 = @(x1) (exp(-log(x2(x1) + W12(2).*x1(i)) - (x1(i).*(W12(1)./(x1(i) + W12(1).*x2(x1)) - W12(2)./(W12(2).*x1(i)+x2(x1))))));
    
    P1sat = @(T) 10.^(A1(1) - A1(2)./(T + A1(3)));
    P2sat = @(T) 10.^(A2(1) - A2(2)./(T + A2(3)));
    
    y1 = @(T) gamma1(x1).*(P1sat(T).*x1(i)./P);
    y2 = @(T) gamma2(x1).*(P2sat(T).*x2(x1)./P);
    
    totY = @(T) y1(T) + y2(T) - 1;
    
    Tb(i) = fsolve(totY, (Tb1+Tb2)/2, options);
    
    Y(i) = y1(Tb(i));
end

%Plotting compositions against different things
plot(x1, Tb)
hold on
plot(Y, Tb)
hold on
axis([0 1 min(Tb) max(Tb)])
xlabel('x1, y1')
ylabel('T [K]')
grid on

figure(2)
plot(x1, Y)
hold on
plot(x1, x1)
xlabel('x1')
ylabel('y1')
axis([0 1 0 1])
hold off

%Molar flows
F = 34.099; D = 10; B = F - D;
xF = 0.77; xD = 0.6; xB = (xF*F-xD*D)/B;

q = 1;    %Liquid, ignore this value

yF = nonidealTb(P, Tb1, A1, A2, xF, W12);
LV = (xD - yF)/(xD - xF);
Rmin = LV/(1 - LV)
R = 1.5*Rmin

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V;
Lbar = L + F;

%Initial temperature estimation and starting composition
Tstart = (Tb1 + Tb2)/2;
x(1) = xB;
y(1) = nonidealTb(P, Tstart, A1, A2, x(1), W12);

%%
%Infinite loop in nonidealTb because of Antoine constants? Works with other groups values Eva notes:
%it it shitty poo poo

%Stepping through the stripping section of the distillation tower
i = 1; 
while x<xF
    
    i = i + 1;
    x(i) = Vbar/Lbar*y(i-1) + B/Lbar*xB;
    y(i) = nonidealTb(P, Tstart, A1, A2, x(i), W12);
    
end

%%

%Stepping through the enriching section of the distillation tower
while y<xD
    
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F*xF);
    y(i) = nonidealTb(P, Tstart, A1, A2, x(i), W12);
    i = i + 1;
    
end

idealTrays = length(x) - 1    %Ideal tray's ignoring the reboiler
realTrays = ceil(idealTrays/0.7)    %Actual amount of trays given a 70% efficiency and rounded up