
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
% differential equation solver where ode_func is the differential function
% at hand and the rest is the parameters and cat & Y matrixes. 

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(cat,XA,''), hold on               % plot with conversion against catalyst mass.
xlabel('Amount of catalyst (kg)')
ylabel('XA (%)')

figure(2)
set ( gca, 'xdir', 'reverse' )            %x axis get set to reverse number order (largest to smallest)
plot(T,XA,''), hold on                   %plot with conversion against temperature 
xlabel('Temperature (K)')
ylabel('XA (%)')

Tslut(e+1)=T(end);                         %end temperatures of each reactor for future use
XA_start=max(XA);  %set the new start conversion for the new reactor as the end conversion of the reactor in this itteration
e=e+1;                                                                     %    counting the number of reactors
leg(e,:)= "Reaktor "+ e;                                           %    matrix for the legend of the graphs
end
figure(1)
legend(leg,'location','southeast') %legend with specific placement
figure(2)
legend(leg,'location','northeast') % -----------  I  I  ---------------
disp("reaktorer: "+e)                    % displaying the number of reactors.






