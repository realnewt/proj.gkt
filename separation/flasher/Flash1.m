%Antoine constants
clear
clc
clf
A1 = [15.5381, 2032.73, -33.15]; %butan
A2 = [18.3036, 3816.44, -46.13]; %vatten

%total pressure
P = 760;  %mmHg

x1 = linspace(0, 1);
options = optimset('Display', 'off');
%For loop which will calculate the temperature which gives y1+y2=1 for each x1 for the
%system based on Raoults law and assumption of ideal mixture
for i=1:length(x1)
x2=@(x1)(1-x1(i));

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3)))); %calculating Psat for (1) and (2)

y1=@(T)(P1sat(T).*x1(i)./P); y2=@(T)(P2sat(T).*x2(x1)./P); %using Raoults law to calculate y1 and y2

ysum=@(T)(y1(T)+y2(T)-1); %creating the function to be solved 

T(i)=fsolve(ysum,300,options); %finding the T which gives ysum=0

y11(i)=y1(T(i)); %calculating these ys and saving in a vector using the functions created in line 18 and 20
end

%plotting the temperature vs x-y for the system
plot(x1,T,'blue')
hold on
plot(y11,T,'red')
axis([0 1 min(T) max(T)])
xlabel('Fraction of "Butane"')
ylabel('T [K]')
title('Temperature vs composition for water/"butane" system')
grid on
legend('Liquid phase composition', 'Vapor phase composition')

%% Flash 1
T=350; %temperature at which the flash will operate
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T); %function to be solved
x11=fzero(fun,0.1) %solving the function to find the component fraction of the liquid phase
[T,y1]=find_Tb(P,T,A1,A2,x11) %using the found x-value to calculate the component fraction in the vapour phase

Fin = [11.5002, 42.4998, 42.4998, 540]; %Molflöden: Buthane|Buthene|Hydrogen|Water
Ftot=sum(Fin); 
z1in=sum(Fin(1:3))/Ftot; %defining the inlet component fraction 
A=[1 1 Ftot;y1 x11 Ftot*z1in]; %assembling the equations to be solved
flows=rref(A); %solving the equations by row reduction. First line gives flow for vapour, second gives liquid flow
line([z1in z1in],[0 T]); line([x11 y1], [T T]); line([x11 x11], [0 T], 'Linestyle', '--'); line([y1 y1], [0 T], 'Linestyle', '--') %drawing the flash operation into the diagram
legend('Liquid phase composition', 'Vapor phase composition')
%% Flash 2
%repeating the same for flash 1 but now at a new operating temperature
T2=320;
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T2);
P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3))));
x12=fzero(fun,0.2)
[T,y12]=find_Tb(P,T,A1,A2,x12)

Ftot=flows(1,3);
A2=[1 1 Ftot; y12 x12 y1*Ftot]
flows2=rref(A2) 

line([y1 y1],[0 T2]); line([x12 y12], [T2 T2]); line([x12 x12], [0 T2], 'Linestyle', '--'); line([y12 y12], [0 T2], 'Linestyle', '--')
%% Flash 3
T3=290;
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T3);
P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3))));
x13=fzero(fun,0.4)
[T,y13]=find_Tb(P,T,A1,A2,x13)

Ftot=flows2(1,3);
A3=[1 1 Ftot; y13 x13 y12*Ftot]
flows3=rref(A3) 

line([y12 y12],[0 T3]); line([x13 y13], [T3 T3]); line([x13 x13], [0 T3], 'Linestyle', '--'); line([y13 y13], [0 T3], 'Linestyle', '--')
legend('Liquid phase composition', 'Vapor phase composition')
clc;
%% Flow after first flash
butout1=y1*flows(1,3); %the flow out of the butane component in the vapour phase
Flow_1=[Fin(1)/sum(Fin(1:3))*butout1; %calculating components to be the new feed into the second flash
        Fin(2)/sum(Fin(1:3))*butout1;
        Fin(3)/sum(Fin(1:3))*butout1;
        (1-y1)*flows(1,3)]

%% Flow after second flash
butout2=y12*flows2(1,3); 


Flow_2=[Fin(1)/sum(Fin(1:3))*butout2;
        Fin(2)/sum(Fin(1:3))*butout2;
        Fin(3)/sum(Fin(1:3))*butout2;
        (1-y12)*flows(1,3)]

%% Final flows out
butout3=y13*flows3(1,3); 


Flow_3=[Fin(1)/sum(Fin(1:3))*butout3; %final components out
        Fin(2)/sum(Fin(1:3))*butout3;
        Fin(3)/sum(Fin(1:3))*butout3;
        (1-y13)*flows(1,3)]

