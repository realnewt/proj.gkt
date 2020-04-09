% Program to solve a system of coupled ordinary differential equations
% 
% Reaction carried out in Batch Reactor
% Reaction:        A + 2B -> P
% Rate Epression:  r = k * CA * CB^2 / (1 + K * CP)
% 
clear 
clc
%Set values of constants
FA0=10;  FB0=50; T0=850+273;  % initial number of moles of A, B and P
R=8.314;    T=T0;       Ea=141e3;
k=0.0596;  K1=22.90; K2=7.56 ; Ke=2.1e7*exp(-122 /(R*T));           % kinetic par. k @550 K1&K2 bar Ke scale with temp also kJ/mol bar
                     % volume of reactor
V=0
g=0;
T=550+273.15;   A=k/(exp(-Ea/(R*T)));  


%Start and Stop times
V_start=0;   % start value for t
V_final=100; % final value for t
%[t_out,U_out] = ode45(@ode_function,[t_range],[initial_values], [options],constants)
[V, U]=ode15s(@ode_eq,[V_start, V_final],[FA0 FB0 T0],[],k,Ke,K1,K2,V,A,Ea,R);
%extract NA, NB and NP from U
FA=U(:,1);  FB=U(:,2); 
%calculate concentrations


%calculate conversion of A
XA=(FA0-FA)./FA0

%% plot moles versus time
figure(1)
plot(t,NA,t,NB,t,NP)
xlabel('t')
ylabel('moles')
%plot conversion of A versus time
figure(2)
plot(t,XA)
xlabel('t')
ylabel('Conversion of A')