function dYdW=ode_func(cat,Y,HR,P,CP,FA0,FB0,FW0)
%function file containing differential equations

XA=Y(1); T=Y(2);

R=8.314; %J/K/mol
K1=22.90; %bar^-1.5
K2=7.56; %/bar
Ea=141e3; %J/mol

FA=(FA0)*(1-XA); FB=FB0+FA0*XA; FH=FA0*XA; FW=FW0;
Ftot=FA+FB+FH+FW;

PA=FA/Ftot*P; PB=FB/Ftot*P; PH=FH/Ftot*P;
k=0.0596*exp((Ea/R).*(1./(550+273)-1./T)); %mol/(kg cat)/s/bar
Ke=2.1*10^7*exp(-122/(R*T)); %bar

cpA=@(T)CP(1,1)+CP(1,2).*T+CP(1,3).*T.^2+CP(1,4).*T.^3;
cpB=@(T)CP(2,1)+CP(2,2).*T+CP(2,3).*T.^2+CP(2,4).*T.^3;
cpH=@(T)CP(3,1)+CP(3,2).*T+CP(3,3).*T.^2+CP(3,4).*T.^3;
cpW=@(T)CP(4,1)+CP(4,2).*T+CP(4,3).*T.^2+CP(4,4).*T.^3;

delHr=HR-integral(cpA,298,T)+integral(cpB,298,T)+integral(cpH,298,T);

r=(k*(PA-(PB*PH)/Ke))/(1+K1*PB*PH^0.5+(K2*PH)^0.5);

dYdW=[r/FA0 %mole balance
    (r*(-delHr))/(FA*cpA(T)+FB*cpB(T)+FH*cpH(T)+FW*cpW(T))]; %heat balance