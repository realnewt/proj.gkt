clc,clear
M_CO=0.02801;      M_H2=0.00201568;    M_H2O=0.01801528;   M_CH3OH=0.03104;

a_CO=30.87;      b_CO=-0.01357;     c_CO=2.789*10^-5;       d_CO=-1.272*10^-8;
a_H2=27.14;      b_H2=0.009274;     c_H2=-1.381*10^-5;      d_H2=7.645*10^-9;
av_H2O=72.43;    bv_H2O=0.01039;    cv_H2O=-1.497*10^-6;    dv_H2O=0;  % index v står för vätska
ag_H2O=32.24;    bg_H2O=0.001924;   cg_H2O=1.055*10^-5;     dg_H2O=-3.596*10^-9;   % index  g står för gas
av_CH3OH=11.7;   bv_CH3OH=-0.4264;  cv_CH3OH=0.00109;       dv_CH3OH=0;
ag_CH3OH=21.15;  bg_CH3OH=0.07092;  cg_CH3OH=2.587*10^-5;   dg_CH3OH=-2.852*10^-8;

%  molflöde  + massa värme växlare 1
F1_CO=5398.3;             F1_H2=10797;              F1_H20=1799.4;          
m1_CO=(F1_CO*M_CO)/60;    m1_H2=(F1_H2* M_H2)/60;   m1_H2O=(F1_H20* M_H2O)/60;  
T0=283;         T0a=400;     T1_medel=(T0+T0a)/2;
Cp1_CO=Cp(T1_medel,a_CO,b_CO,c_CO,d_CO,M_CO);
Cp1_H2=Cp(T1_medel,a_H2,b_H2,c_H2,d_H2,M_H2);

q1_CO=m1_CO*Cp1_CO*(T0a-T0);
q1_H2=m1_H2*Cp1_H2*(T0a-T0);
q1_gas=q1_CO+q1_H2

%  vatten genomgår fasomvandling
T0=283;         T01=373;     T12_medel=(T0+T01)/2;
Cp1v_H2O=Cp(T12_medel,av_H2O,bv_H2O,cv_H2O,dv_H2O,M_H2O);
q1_H2O_liquid=m1_H2O*Cp1v_H2O*(T01-T0);

T01=373;         T0a=300+273;     T13_medel=(T01+T0a)/2;
Cp1g_H2O=Cp(T13_medel,ag_H2O,bg_H2O,cg_H2O,dg_H2O,M_H2O);
q1_H2O_gas=m1_H2O*Cp1g_H2O*(T0a-T01);
q1_H2O=q1_H2O_liquid+q1_H2O_gas   %+2258000   Hvap gör att svaret blir orimlig

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


