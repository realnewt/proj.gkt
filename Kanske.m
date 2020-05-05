
% A=butane B=butene H=hydrogen W=water
% matlabprogram to determine the number of reactors in the reaction 
%Butane --> Butene+H2 with assumed amount of katalyst and heating in
%between the reactors.
clear, clc;
%clf
format short
e=0;
cat_tot=500; %kg cat
XA_start=0;
T0=950; %constant start temp. for every reactor in kelvin
FA0=54; FB0=0.5; FW0=10*FA0; FH0=0; %molar flowrates in to the first reactor mol/s

HR=116.3e3; %J/mol reaction enthalpy
P=1; %bar constant total pressure
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09;
    32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
    
Tslut=zeros(1,5); %define a matrix just so matlab doesn't complain
XAut=zeros(1,5);
while XA_start-0.78<1e-4
%{
a while-loop to determine the amount of reactors needed 
to get the conversion sought after. the loop stops when converison hits a
specific level (0.95 in this case).
%}
[cat,Y]=ode15s(@(cat,Y) ode_func(cat, Y, HR, P, CP, FA0,FH0, FB0, FW0), [0 cat_tot], [XA_start T0]); 
% differential equation solver where ode_func is the differential function
% at hand and the rest is the parameters and cat & Y matrixes. 

XA=Y(:,1); T=Y(:,2);
%------------------------------
%{
figure(1)
plot(cat,XA,'linewidth',1), hold on               % plot with conversion against catalyst mass.
xlabel('Mängd katalysator (kg)')
ylabel('XA')

figure(2)
set ( gca, 'xdir', 'reverse' )            %x axis get set to reverse number order (largest to smallest)
plot(T,XA,'linewidth',1), hold on                   %plot with conversion against temperature 
xlabel('Temperatur (K)')
ylabel('XA')
%} 
%-----------------------------------------

Tslut(e+1)=T(end); %end temperatures of each reactor for future use 
XA_start=max(XA);  %set the new start conversion for the new reactor as the end conversion of the reactor in this itteration
fa=(FA0)*(1-XA_start);
fb=FB0+FA0*XA_start; fh=FH0+FA0*XA_start; fw=FW0;
ftot=fa+fb+fh+fw;
XAut(e+1)=XA_start(end); 
e=e+1;                  %    counting the number of reactors
ff(e,:)=[fa fb fh fw];
leg(e,:)= "Reaktor "+ e;                                           %    matrix for the legend of the graphs
end
figure(1)
legend(leg,'location','southeast') %legend with specific placement
figure(2)
legend(leg,'location','northeast') % -----------  I  I  ---------------
disp("reaktorer: "+e)                    % displaying the number of reactors.

disp('grattis')


%-----------------------------------------------------



% Molmassa för komponenterna
M=[58.12*10^-3 ;
   56.1*10^-3 ;
   2*1.00784*10^-3 ;
   18.01528*10^-3];                             %Butan-Buten-H2-H2O i kg/mol


for i=1:4
    %  molflöde  + massa värme växlare 1
F_mol=[ff(1,i);
       ff(2,i);
       ff(3,i);
       ff(4,i)];
   F_massa=[ff(1,i)*M(1);
         ff(2,i)*M(2);
         ff(3,i)*M(3);
         ff(4,i)*M(4);];                    %Butan-Buten-H2-H2O i kg/s  
     T0=Tslut(i); T0a=950;   T=(T0+T0a)/2;
cpBA=@(T)CP(1,1)+CP(1,2).*T+CP(1,3).*T.^2+CP(1,4).*T.^3; cpBA=cpBA(T);
cpBE=@(T)CP(2,1)+CP(2,2).*T+CP(2,3).*T.^2+CP(2,4).*T.^3; cpBE=cpBE(T);
cpH=@(T)CP(3,1)+CP(3,2).*T+CP(3,3).*T.^2+CP(3,4).*T.^3; cpH=cpH(T);
cpW=@(T)CP(4,1)+CP(4,2).*T+CP(4,3).*T.^2+CP(4,4).*T.^3; cpW=cpW(T);

QBA=F_massa(1)*cpBA*(T0a-T0);
QBE=F_massa(2)*cpBE*(T0a-T0);
QH=F_massa(3)*cpH*(T0a-T0);
QW=F_massa(4)*cpW*(T0a-T0);
end     


%Cp=Cp(T1_medel,a_CO,b_CO,c_CO,d_CO,M_CO);
%Cp1_H2=Cp(T1_medel,a_H2,b_H2,c_H2,d_H2,M_H2);


%q1_CO=m1_CO*Cp1_CO*(T0a-T0);
%q1_H2=m1_H2*Cp1_H2*(T0a-T0);
%q1_gas=q1_CO+q1_H2

%  vatten genomgår fasomvandling
T0=283;         T01=373;     T12_medel=(T0+T01)/2;
Cp1v_H2O=Cp(T12_medel,av_H2O,bv_H2O,cv_H2O,dv_H2O,M_H2O);
q1_H2O_liquid=m1_H2O*Cp1v_H2O*(T01-T0);

T01=373;         T0a=300+273;     T13_medel=(T01+T0a)/2;
Cp1g_H2O=Cp(T13_medel,ag_H2O,bg_H2O,cg_H2O,dg_H2O,M_H2O);
q1_H2O_gas=m1_H2O*Cp1g_H2O*(T0a-T01);
q1_H2O=q1_H2O_liquid+q1_H2O_ga;s   %+2258000   Hvap gör att svaret blir orimlig

U_gas=50;     C=0.8;    U_liquid=885;   %  Givna  

% inlopps och utlopps tempratur för värmevatten för värmeväxlare 1
T1_vapor1=210+273;
T1_vapor2= T1_vapor1 - C*(T0a-T0);              % C=0.8=(Th1-Th2)/(Tc2-Tc1)
T1_vapor_medel= (T1_vapor1+T1_vapor2)/2;
Cp1_vapor=Cp(T1_vapor_medel,ag_H2O,bg_H2O,cg_H2O,dg_H2O,M_H2O);  
 

%                        &&   VÄRMEVÄXLARE  1    &&

% ångflöde samt area som krävs för upphettning av CO
C1_gas=m1_CO*Cp1_CO+  m1_H2*Cp1_H2;                    % C1=Cmin
Cmax1_gas=C1_gas/0.8;               %                  C=0.8=Cmin/Cmax  

C1_liquid=m1_H2O*Cp1g_H2O+m1_H2O*Cp1v_H2O;
q1=q1_H2O+q1_gas;
C1=C1_gas+C1_liquid;
%Epsilon1=(q1_CO)/((m1_CO*Cp1_CO)*(T1_vapor1-T0));  %      epsilon=q/(Cmin(Th_in-Tc_in)
Epsilon1=q1/(C1*(T1_vapor1-T0));  %      epsilon=q/(Cmin(Th_in-Tc_in)


Cmax1_liquid=C1_liquid/0.8;               



m1_vapor_gas=Cmax1_gas/Cp1_vapor;
m1_vapor_liquid=Cmax1_liquid/Cp1_vapor;
m1_vapor=m1_vapor_gas+ m1_vapor_liquid       % required heat-medium



tdrift=8000;   
beta=0.16*10^-3;              %[SEK/Wh]       Kostnad för ångan
A1_gaser=Area1(Epsilon1,C,C1,U_gas) 
A1_liquid=Area1(Epsilon1,C,C1,U_liquid) 
A1=   A1_gaser + A1_liquid


%%

F3_CO=12596;            F3_H2=25192;    
m3_CO=F3_CO*M_CO;       m3_H2=F3_H2* M_H2;    
T4=274;    T5=523;      T3_medel=(T4+T5)/2;
Cp3_CO=Cp(T3_medel,a_CO,b_CO,c_CO,d_CO,M_CO);
Cp3_H2=Cp(T3_medel,a_H2,b_H2,c_H2,d_H2,M_H2);

% inlopps och utlopps tempratur för värmevatten för värmeväxlare 3
T3_vapor1=210+273;     T3_vapor2= T3_vapor1 - C*(T5-T4);         
T3_vapor_medmel= (T3_vapor1+T3_vapor2)/2;
Cp3_vapor=Cp(T3_vapor_medmel,ag_H2O,bg_H2O,cg_H2O,dg_H2O,M_H2O);  

q3_CO=m3_CO*Cp3_CO*(T5-T4);
q3_H2=m3_H2*Cp3_H2*(T5-T4);


%                        &&  VÄRMEVÄXLARE 3    &&
U3=50;
% ångflöde samt area som krävs för upphettning av CO
q3=q3_CO+q3_H2;
C3=m3_CO*Cp3_CO + m3_H2*Cp3_H2;    Cmax3=C3/0.8;    
Epsilon3=q3/(Cmax3*(T3_vapor1-T4));     
A3_vaxlare=Area1(Epsilon3,C,C3,U3)    
m3=Cmax3/Cp3_vapor 

%Ka3=32000+70*(A3_vaxlare)^1.2;   %[SEK/(m2 år)]
%A3_optimal=optimalarea(100,C3,Cmax3,Ka3,T3_vapor1,T4,tdrift,beta,100);
%%


R=8.314;              R_CO=R/M_CO;           R_H2=R/M_H2;   R_H2O=R/M_H2O;
Cv1_CO=Cp1_CO-R_CO;      kappa1_CO=Cp1_CO/Cv1_CO;
Cv1_H2=Cp1_H2-R_H2;      kappa1_H2=Cp1_H2/Cv1_H2;  
Cv1_H2O=Cp1g_H2O-R_H2O;      Kappa1_H2O=Cp1g_H2O/Cv1_H2O;  %Vatten är vätska => pump
kappa1=(kappa1_CO+kappa1_H2+Kappa1_H2O)/3;
P_in=101325;
P_ut=101325*100;
Ctot=m1_CO*Cp1_CO+m1_H2*Cp1_H2+ m1_H2O*Cp1g_H2O;   % C1=Cmin

komp_Area=kompressor(Ctot,kappa1,P_in,T0a,P_ut,0.8)

Tin=400;
eta_is=0.8;
P_step = (P_ut/P_in)^(1/3);  %[]
%Temperatur ut från varje kompressorsteg för isentrop kompression.
Tut_is = Tin*P_step^((kappa1-1)/kappa1);  %[K] 
%Verklig temperatur ut från varje kompressorsteg.
Tut = Tin + (Tut_is-Tin)/eta_is %[K] 
