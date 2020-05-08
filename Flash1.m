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
for i=1:length(x1) 
x2=@(x1)(1-x1(i));

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3))));

y1=@(T)(P1sat(T).*x1(i)./P); y2=@(T)(P2sat(T).*x2(x1)./P);

ysum=@(T)(y1(T)+y2(T)-1);

T(i)=fsolve(ysum,300,options);

y11(i)=y1(T(i));
end

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
T=350;
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T);
P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3))));
x11=fzero(fun,0.1)
[T,y1]=find_Tb(P,T,A1,A2,x11)

Fin = [11.5002, 42.4998, 42.4998, 540]; %Molflöden: Buthane|Buthene|Hydrogen|Water
Ftot=sum(Fin); 
z1in=sum(Fin(1:3))/Ftot;
A=[1 1 Ftot;y1 x11 Ftot*z1in];
flows=rref(A) %utflöden från den första flashen, vapor översta raden
line([z1in z1in],[0 T]); line([x11 y1], [T T]); line([x11 x11], [0 T], 'Linestyle', '--'); line([y1 y1], [0 T], 'Linestyle', '--')
legend('Liquid phase composition', 'Vapor phase composition')
%% Flash 2
T2=320;
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T2);
P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3))));
x12=fzero(fun,0.2)
[T,y12]=find_Tb(P,T,A1,A2,x12)

Ftot=flows(1,3);
A2=[1 1 Ftot; y12 x12 y1*Ftot]
flows2=rref(A2) %utflöden från den andra flashen, vapor först

line([y1 y1],[0 T2]); line([x12 y12], [T2 T2]); line([x12 x12], [0 T2], 'Linestyle', '--'); line([y12 y12], [0 T2], 'Linestyle', '--')
%% Flash 3
T3=290;
fun=@(x1)(find_Tb(P,T,A1,A2,x1)-T3);
P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3))));
x13=fzero(fun,0.4)
[T,y13]=find_Tb(P,T,A1,A2,x13)

Ftot=flows2(1,3);
A3=[1 1 Ftot; y13 x13 y12*Ftot]
flows3=rref(A3) %utflöden från den tredje flashen vapor först

line([y12 y12],[0 T3]); line([x13 y13], [T3 T3]); line([x13 x13], [0 T3], 'Linestyle', '--'); line([y13 y13], [0 T3], 'Linestyle', '--')
legend('Liquid phase composition', 'Vapor phase composition')

butout3=y13*flows3(1,3) %slutgiltiga utflödet av butan som är att jämför med
butin1=sum(Fin(1:3)) %inflödet av butan till första flashen