clear, clc, clf

%1 = buten, 2 = butan

%Wilson parameters
W12 = 0.48584; 
W21 = 1.64637;

%Antoine constants for degC, mmHg, log10

A2=15.7564; B2=2132.42; C2=-33.15 ;%buten
A1=15.5381; B1=2032.73; C1=-33.15;%butan
%total pressure
P =3.3*760;  %mmHg

tb1=B1/(A1-log(P))-C1;
tb2=B2/(A2-log(P))-C2;

x1 = linspace(0,1,1000);
Tstart=(tb1+tb2)/2;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);

    %activity coefficients at x1
    gamma2 = exp(-log(x1(i)+W12*x2)+x2.*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
    gamma1 = exp(-log(x2+W21*x1(i))-x1(i).*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
   
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tb(T,x1(i),x2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = exp(A1-B1./(Tb(i)+C1));
    P1 = gamma1.*P01.*x1(i);
    y1(i) = P1./P;
end
figure(1)
hold on
axis([0 1 min(Tb) max(Tb)])
plot(x1, Tb)
plot(y1, Tb)
xlabel('x1,y1')
ylabel('T [K]')

hold off

figure(2)
hold on
plot(x1,y1)
plot(x1,x1,'red')
xlabel('x')
ylabel('y')
%%
%Molar flows
%F = 34.099; 
F = 34.099; xF=0.23; xD=0.45; xB=0.1; %Känt inflöde F o xF, önskade sammansättningar xD(toppen) och xB=azeotrop(botten)
A=[1 1 F; xD xB xF*F]; Flows=rref(A); %räknar ut toppflöde och bottenflöde med total samt komponentbalans
D=Flows(1,3); B=Flows(2,3);
q = 1;

%Calculate yF
Tbf = fsolve(@(T)find_Tb(T,xF,1-xF,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tbf+C1));
P1 = gamma1.*P01.*xF;
yF= P1./P; 

%Calculate R
LV = (xD - yF)/(xD - xF);
Rmin = LV/(1 - LV);
R = 2*Rmin;

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V;
Lbar = L + F;

%Initial temperature estimation and starting composition bottom of tower
Tstart = (tb1 + tb2)/2;
x(1) = xB;
Tb = fsolve(@(T)find_Tb(T,xB,1-xF,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tb+C1));
P1 = gamma1.*P01.*xB;
y(1)= P1./P;
%% Botten
i = 1; 
while x<xF
    i = i + 1;
    x(i)=Vbar/Lbar*y(i-1) + B/Lbar*xB;
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21)
    
end
%% toppen
while y<xD
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F*xF);
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21)
    i = i + 1;
end
n=i-1;
real=n/0.7
%% Torndimensioner

ts=0.45;
H=(real+1)*ts

%Medelmolmassor
ML=x(1)*58.12+(1-x(1))*56.11;
MV=y(1)*58.12+(1-y(1))*56.11;

%Temperatur i botten
Tbb = fsolve(@(T)find_Tb(T,x(1),1-x(1),gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options)

%Medeldensiteter

MVrho=((ML/1000)*(P*133.322368))/(8.314*Tbb); %Kg/m^3
MLrho=(x(1))*563+(1-x(1))*588; %Kg/m^3

FLV=(Lbar*ML)/(Vbar*MV)*sqrt(MVrho/MLrho);
CF=0.28; %flooding constant from diagram
sigma=15.8; %taken from some sketchy page isobutene at 20Celcius
FST=(sigma/20)^0.2;
C=CF*FST;
Uf=C*sqrt(((MLrho-MVrho)/MVrho));
ada=0.1+(FLV-0.1)/9;
DT=sqrt((4*V*(MV/1000))/(0.8*(Uf/3.28)*pi*(1-ada)*MVrho))

%% Totalkondensor och återkokare
Hvap1=21.5*1000 %J/mol isobutan
Hvap2=22.5*1000 %J/mol isobutene
Havgtop=xD*Hvap1+(1-xD)*Hvap2
Havgbot=x(1)*Hvap1+(1-x(1))*Hvap2

%condenser

Qc=D*(R+1)*Havgtop %Joule/s

%reboiler
Qr=Vbar*Havgbot %Joule/s










