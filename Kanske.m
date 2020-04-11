% A=butan B=buten H=väte W=vatten
clear
clc
clf
Wtot=20000; %kg cat
XA_start=0;
T0=1000; %K
FA0=100; %mol/s
FB0=0;
FW0=10*FA0;

HR=116.3e3; %J/mol
P=1; %bar
CP=[1.39 0.3847 -1.846e-04 2.895e-08; 16.05 0.2804 -1.091e-04 9.098e-09; 27.14 0.009274 -1.3813e-05 7.645e-09; 32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0),[0 Wtot],[XA_start T0]);

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(W,XA, 'blue'), hold on
xlabel('amount catalyst(kg)')
ylabel('XA')

figure(2)
plot(T,XA, 'blue'), hold on
xlabel('Temp(K)')
ylabel('XA')


%%

Wtot=20000; %kg cat
XA_start=0;
T0=1000; %K
FA0=100; %mol/s
FB0=.5*FA0;
FW0=10*FA0;

HR=116.3e3; %J/mol
P=1; %bar
CP=[1.39 0.3847 -1.846e-04 2.895e-08; 16.05 0.2804 -1.091e-04 9.098e-09; 27.14 0.009274 -1.3813e-05 7.645e-09; 32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0),[0 Wtot],[XA_start T0]);

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(W,XA, 'r')
xlabel('amount catalyst(kg)')
ylabel('XA')

figure(2)
plot(T,XA, 'r')
xlabel('Temp(K)')
ylabel('XA')






