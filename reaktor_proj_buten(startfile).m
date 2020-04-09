% Program to solve a system of coupled ordinary differential equations
% 
% Reaction carried out in Batch Reactor
% Reaction:        A + 2B -> P
% Rate Epression:  r = k * CA * CB^2 / (1 + K * CP)
% 
clear all;
%Set values of constants
NA0=10;  NB0=50;  NP0=0; T0=175+273;  % initial number of moles of A, B and P
k=0.005;  K=0.5;           % kinetic parameters
V=10;                     % volume of reactor
%Start and Stop times
t_start=0;   % start value for t
t_final=100; % final value for t
%[t_out,U_out] = ode45(@ode_function,[t_range],[initial_values],
[options],constants)
[t U]=ode15s(@ode_eq,[t_start t_final],[NA0 NB0 NP0 T0],[],k,K,V);
%extract NA, NB and NP from U
NA=U(:,1);  NB=U(:,2);  NP=U(:,3);
%calculate concentrations
CA=NA/V;  CB=NB/V;  CP=NP/V;
%calculate conversion of A
XA=(NA0-NA)./NA0;
%plot moles versus time
figure(1)
plot(t,NA,t,NB,t,NP)
xlabel('t')
ylabel('moles')
%plot conversion of A versus time
figure(2)
plot(t,XA)
xlabel('t')
ylabel('Conversion of A')