function dUdW=ode_eq(W,U,k,Ke,K1,K2,V,A,Ea,R)
% function file containing differential equations
%extract NA, NB and NP from U
FA=U(1);  FB=U(2);  T=U(3);

%calculate concentrations preassure
Pa=0;    Pb=0;      PH=0;
%calculate reaction rate
k=A*exp(-Ea/R/T);
r=(k*(Pa-(Pb*PH/Ke)))/(1+K1*Pb*PH^.5+(K2*PH)^.5);
%differential equations
dFAdw=-r;
dFBdW=r;
dTdW=0;      %  (r*(-dH)+Q)/(Fa*Cpa+Fb*Cpb);
%Assemble Differential Equations
dUdW=[dFAdw
    dFBdW
    dTdW];