function dUdW=ode_eq(W,U,k,K,V)
% function file containing differential equations
%extract NA, NB and NP from U
FA=U(1);  FB=U(2);  T=U(3);

%calculate concentrations preassure
Pa=;    Pb=;      PH=;
%calculate reaction rate
k=A*exp(-Ea/R/T);
r=(k(Pa-(PbPH/Ke)))/(1+K1*Pb*PH^.5+(K2PH)^.5);
%differential equations
dFAdw=-r;
dFBdW=r;
dTdW=0;
%Assemble Differential Equations
dUdW=[dFAdw
    dFBdW
    dTdW];
