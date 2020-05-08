%Antoine constants
clear
clc
clf

A1 = [13.6333, 164.90, 3.19]; %hydorgen
A2 = [15.5381, 2032.73, -33.15]; %butan

%total pressure
P = 4*760;  %mmHg

x1 = linspace(0, 1);
options = optimset('Display', 'off');
%For loop which will calculate the temperature which gives y1+y2=1 for each x1 for the
%system based on Raoults law and assumption of ideal mixture
for i=1:length(x1) 
x2=@(x1)(1-x1(i));

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3)))); %calculating Psat for (1) and (2)

y1=@(T)(P1sat(T).*x1(i)./P); y2=@(T)(P2sat(T).*x2(x1)./P); %using Raoults law to calculate y1 and y2

ysum=@(T)(y1(T)+y2(T)-1); %creating the function to be solved 

T(i)=fsolve(ysum,200,options);%finding the T which gives ysum=0

y11(i)=y1(T(i)); %calculating these ys and saving in a vector using the functions created in line 18 and 20
end

%plotting the temperature vs x-y for the system
plot(x1,T,'blue')
hold on
plot(y11,T,'red')
axis([0 1 min(T) max(T)])
xlabel('Fraction of "Hydrogen"')
ylabel('T [K]')
title('Temperature vs composition for hydrogen/"butane" system')
grid on
legend('Liquid phase composition', 'Vapor phase composition')
%%
T=155; %temperature at which the flash will operate
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T);  %function to be solved
x11=fzero(fun,0.1) %solving the function to find the component fraction of the liquid phase
[T,y1]=find_Tb(P,T,A1,A2,x11) %using the found x-value to calculate the component fraction in the vapour phase

Fin = [7.1328, 26.36, 26.36, 0.25591]; %molar flows: Buthane|Buthene|Hydrogen|Water
Ftot=sum(Fin); 
z1in=Fin(3)/Ftot; %defining the inlet component fraction 
A=[1 1 Ftot;y1 x11 Ftot*z1in];%assembling the equations to be solved
flows=rref(A); %solving the equations by row reduction. First line gives flow for vapour, second gives liquid flow
line([z1in z1in],[0 T]); line([x11 y1], [T T]); line([x11 x11], [0 T], 'Linestyle', '--'); line([y1 y1], [0 T], 'Linestyle', '--') %drawing the flash operation into the diagram
legend('Liquid phase composition', 'Vapor phase composition')

%calculating components flow out of the reactor to go to the distillation
%tower
butaneout=Fin(1)/(Fin(1)+Fin(2)+Fin(4))*(1-x11)*flows(2,3)
buteneout=Fin(2)/(Fin(1)+Fin(2)+Fin(4))*(1-x11)*flows(2,3)
waterout=Fin(4)/(Fin(1)+Fin(2)+Fin(4))*(1-x11)*flows(2,3)
hydrogenout=x11*flows(2,3)



