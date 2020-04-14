
% A=butane B=butene H=hydrogen W=water
% matlabprogram to determine the number of reactors in the reaction 
%Butane --> Butene+H2 with assumed amount of katalyst and heating in
%between the reactors.
clear, clc, clf, format short
e=0;
cat_tot=500; %kg cat
XA_start=0;
T0=950; %constant start temp. for every reactor in kelvin
FA0=54; FB0=0.5; FW0=10*FA0; %molar flowrates in to the first reactor mol/s

HR=116.3e3; %J/mol reaction enthalpy
P=1; %bar constant total pressure
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09;
    32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
    
Tslut=zeros(1,5); %define a matrix just so matlab doesn't complain
while XA_start-0.89<1e-4
%{
a while-loop to determine the amount of reactors needed 
to get the conversion sought after. the loop stops when converison hits a
specific level (0.95 in this case).
%}
[cat,Y]=ode15s(@(cat,Y) ode_func(cat, Y, HR, P, CP, FA0, FB0, FW0), [0 cat_tot], [XA_start T0]);

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(cat,XA, ''), hold on
xlabel('amount catalyst(kg)')
ylabel('XA')

figure(2)
set ( gca, 'xdir', 'reverse' )
plot(T,XA, ''), hold on
xlabel('Temp(K)')
ylabel('XA')

Tslut(e+1)=T(end);

XA_start=max(XA);
FA=(FA0)*(1-XA_start); FB=FB0+FA0*XA_start; 
FH=FA0*XA_start; FW=FW0;
Ftot=FA+FB+FH+FW;
FA0=FA; %mol/s
FB0=FB;
FW0=10*FA0;
e=e+1;
leg(e,:)= "Reaktor "+ e;

end
figure(1)
legend(leg)
figure(2)
legend(leg,'location','northwest')
disp("reaktorer: "+e)






