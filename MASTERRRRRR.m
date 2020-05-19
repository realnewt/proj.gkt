%% Reaktor 1
clear
clc
clf
options = optimset('Display', 'off');
e=0;
Wtot=3000; %kg cat
XA_start=0;
T0=950; %K
FA0=4*54; %mol/s
FB0=0.5;
FH0=0;
FW0=10*FA0;

HR=116.3e3; %J/mol
P=1; %bar
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09; 
    32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0,FH0),[0 Wtot],[XA_start T0]);

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(W,XA)
xlabel('amount catalyst(kg)')
ylabel('XA')

figure(2)
plot(T,XA)
xlabel('Temp(K)')
ylabel('XA')
%kostnad
rhobed=1120; %kg/m3
Volym=Wtot/rhobed;
D=(2*Volym/pi)^(1/3);
L=2*D;

CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
Pi=pi;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)
t=(P*D)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 
V_inner=Pi*(((D/2)-t)^2)*(L-2*t);     %[m^3] Volume of inner tank (air)
Volume_container=Volym-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass_dist.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);  Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B); %[$] Accounts for inflation
Cost_Dist_v_304_single=Cost_Year_A_v_304*10;        %[SEK]
Cost_Dist_h_304_single=Cost_Year_A_h_304*10;        %[SEK]

Cost_Dist_v_304_five=Cost_Dist_v_304_single*1.5;
Cost_Dist_h_304_five=Cost_Dist_h_304_single*1.5;
fprintf('REACTOR ONE\n')
fprintf('Inflow: %2.0f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\n',FA0,FB0,FH0,FW0)
fprintf('Dimensions: diameter is %0.3f m and length is %1.2f m. Volume is %0.2f m^3 with a shell voulme of %0.4f m^3 \n',D,L,Volym,Volume_container)
fprintf('Wall thickness is %0.3f mm and total mass of container is %2.0f kg. Amount of catalyst is %4.0f kg\n',t*1000,Shell_mass_dist,Wtot)
fprintf('Total cost is %f SEK with a final conversion of %0.2f\n',min(Cost_Dist_v_304_five,Cost_Dist_h_304_five),XA(end));
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))

%% Ugn 1
  M=[58.12*10^-3 ;
   56.1*10^-3 ;
   2*1.00784*10^-3 ;
   18.01528*10^-3];

F_massa=[FA0*(1-XA(end))*M(1);
             (FB0+XA(end))*FA0*M(2);
             (FH0+XA(end))*FA0*M(3);
             FW0*M(4);];
         
Tin=T(end); Tut=950; Tmedel=(Tin+Tut)/2;
cpBA=@(Tmedel)CP(1,1)+CP(1,2).*Tmedel+CP(1,3).*Tmedel.^2+CP(1,4).*Tmedel.^3;       cpBA=cpBA(Tmedel);
cpBE=@(Tmedel)CP(2,1)+CP(2,2).*Tmedel+CP(2,3).*Tmedel.^2+CP(2,4).*Tmedel.^3;     cpBE=cpBE(Tmedel);
cpH=@(Tmedel)CP(3,1)+CP(3,2).*Tmedel+CP(3,3).*Tmedel.^2+CP(3,4).*Tmedel.^3;      cpH=cpH(Tmedel);
cpW=@(Tmedel)CP(4,1)+CP(4,2).*Tmedel+CP(4,3).*Tmedel.^2+CP(4,4).*Tmedel.^3;     cpW=cpW(Tmedel);

cp=[cpBA
    cpBE
    cpH
    cpW];

QBA=F_massa(1)*cpBA*(Tut-Tin);
QBE=F_massa(2)*cpBE*(Tut-Tin);
QH=F_massa(3)*cpH*(Tut-Tin);
QW=F_massa(4)*cpW*(Tut-Tin);
Qtot=QBA+QBE+QH+QW;
         
chi=0.8;
Qheat=Qtot/0.8;     %Energi som krävs för uppvärmning J/s
fprintf('OVEN ONE\n')
fprintf('Inflow: %2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n ',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))
fprintf('Energy needed is %f J/s\n', Qheat)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n\n\n ',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,Tut)

%% Reaktor 2
fprintf('REAKTOR 2\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n ',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,Tut)
XA_start=XA(end);
Wtot=3000;
P=1;
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0,FH0),[0 Wtot],[XA_start T0]);
XA=Y(:,1); T=Y(:,2);

figure(3)
plot(W,XA)
xlabel('amount catalyst(kg)')
ylabel('XA')

figure(4)
plot(T,XA)
xlabel('Temp(K)')
ylabel('XA')

rhobed=1120; %kg/m3
Volym=Wtot/rhobed;
D=(2*Volym/pi)^(1/3);
L=2*D;

CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
Pi=pi;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)
t=(P*D)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 
V_inner=Pi*(((D/2)-t)^2)*(L-2*t);     %[m^3] Volume of inner tank (air)
Volume_container=Volym-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass_dist.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);  Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B); %[$] Accounts for inflation
Cost_Dist_v_304_single=Cost_Year_A_v_304*10;        %[SEK]
Cost_Dist_h_304_single=Cost_Year_A_h_304*10;        %[SEK]

Cost_Dist_v_304_five=Cost_Dist_v_304_single*1.5;
Cost_Dist_h_304_five=Cost_Dist_h_304_single*1.5;


fprintf('Dimensions: diameter is %0.3f m and length is %1.2f m\n',D,L)
fprintf('Wall thickness is %0.3f mm and total mass of container is %2.0f kg. Amount of catalyst is %4.0f kg\n',t*1000,Shell_mass_dist,Wtot)
fprintf('Total cost is %f SEK with a final conversion of %0.2f\n',min(Cost_Dist_v_304_five,Cost_Dist_h_304_five),XA(end));
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))

%% Cooler 1  DISCLAIMER: takes a while to run
format shortG
fprintf('COOLER 1\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))
%Data
%Area calculations
C=0.8;      %Temperaturverkningsgrad from kurs PM
U=50;       %[W/m2] KVärmegenomgångstal from kurs PM

%Set up matrix for area calculations
Area_tot=zeros(10,1);

% Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler=[(FA0*(1-XA(end)));
       (FB0+FA0*XA(end));
       (FH0+FA0*XA(end));
       (FW0)];       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass_cooler=[F_mass(F_mol_cooler(1,1),1);
        F_mass(F_mol_cooler(2,1),2);
        F_mass(F_mol_cooler(3,1),3);
        F_mass(F_mol_cooler(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O
A1=15.7564; B1=2132.42; C1=-33.15;%buten
A2=18.3036; B2=3816.44; C2=-46.13; %vatten
% Parameters for adjustment of loop
TH_in=T(end);      %[K] Initial temperature of mixture
TH_ut_final= fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in=273+14;       %Inlet temperature of coolant
A_max=600;      %[m2] Maximum area per heat exchange unit
A_num_max=1000;       %Max number of heat exchange units

%Set up matrices for area calculations
m=zeros(A_num_max,1);       %Set up matrix for TC_out 
p=zeros(A_num_max,1);       %Set up matrix for TH_out
A_per_unit=zeros(A_num_max,1);       %Set up matrix for area per heat exchange unit
r=zeros(A_num_max,1);       %Set up matrix for number of heat exchange units

j=TH_in;      %Adjustment factor for TC_in
k=0;        %Adjustment factor for number of heat exchange units

% Loop calculating number of heat exchange units and area of units
while TH_in>TH_ut_final && k<A_num_max        %While loop that counts number or heat exchange units
k=k+1;
    while j>TH_ut_final       %While loop that iterates TC_in
    j=j-1;        %adjust for finer accuracy 
    
    i=0;
        while i<A_max       %While loop that tests area of each heat exchange unit
        i=i+0.01;        %adjust for finer accuracy    
        
        %Cp calculations mixture
        TH_ut=j;        %Outlet temperature
        TH_medel=(TH_in+TH_ut)/2;       

        CP_matrix=[Cp_new(TH_medel,1,1,1,1,1);
                   Cp_new(TH_medel,2,2,2,2,2);
                   Cp_new(TH_medel,3,3,3,3,3);
                   Cp_new(TH_medel,4,4,4,4,4)];     %[J/Kg*K] Cp each component
     
        %Temperature calculation coolant  
        TC_ut=TC_in+C*(TH_in-TH_ut);        %Outlet temperature of coolant     

        %Energy transfer calculations
        q_matrix=F_mass_cooler(:,1).*CP_matrix(:,1).*(TH_in-TH_ut);      %Energy transfer matrix
        q=sum(q_matrix);        %[J/s] Energy transfer

        %C_max/C_min calculations
        C_max_matrix=CP_matrix.*F_mass_cooler;   
        C_max=sum(C_max_matrix);
        C_min=C*C_max;
        
        %Epsilon calculation
        Epsilon=q/(C_min*(TH_in-TC_in));        %0 < Epsilon < 1
        
        %Function for checking area
        Fun=Epsilon-((1-exp(-(U.*i/C_min)*(1-(C_min/C_max))))/(1-(C_min/C_max)*exp(-(U.*i/C_min)*(1-C_min/C_max))));
        
            if Fun<0        %Checks to see if area is correct
                break;
            end
        end
            if i>A_max      %Checks if max area is exceeded 
                break;
            end
    end
A=i;
T=j;

%Outlet temperature for coolant 
TC_ut=TC_in+(TH_in-TH_ut)./C;

TH_in=T;
j=T;
m(k)=TC_ut;
p(k)=TH_in;
A_per_unit(k)=A;
r(k)=k;
end

% Output of heat exchange loop
Matrix_results=[r,m,p,A_per_unit];        %Matrix with results
column_names={'Number','TC_out','TH_out','Area'};        %Label for matrix with results 
Matrix_results( all(~Matrix_results,2),:)=[];
Heat_exchange=array2table(Matrix_results,'VariableNames',column_names)      %Matrix with results (labeled)

Area_tot=sum(A_per_unit);        %Total area for heat exchange
Area_tot=sprintf('%.f m2',Area_tot);       %Total area of all heat exchange units

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant
Cost_Year_B=zeros(10,1);

Shell_mass=A_per_unit(1,1);      %[m2] Area per unit
S_final=A_per_unit(k,1);        %[m2] Area of final unit

Cost_Year_B=(k-1)*(a_s+b_s*Shell_mass.^n_s)+(a_s+b_s*S_final.^n_s);      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_Cooler=Cost_Year_A*10;      %Conversion from Dollars to SEK

fprintf('Total area needed is %s and total cost is %f SEK \n', Area_tot, Cost_Cooler)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,TH_ut_final)

%% Destillation 1, buten-vatten
clear x
clear y
clear i
clear n
fprintf('DESTILLATION 1\n')
fprintf('Inflow: %2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,TH_ut_final)
%1 = buten, 2 = vatten 

%Antoine constants for degC, mmHg, log10

A1=15.7564; B1=2132.42; C1=-33.15 ;%buten
A2=18.3036; B2=3816.44; C2=-46.13; %vatten
%total pressure
P =1*760;  %mmHg

tb1=B1/(A1-log(P))-C1;
tb2=B2/(A2-log(P))-C2;

x1 = linspace(0,1,1000);
Tstart=(tb1+tb2)/2;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tbnonideal(T,x1(i),x2,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = exp(A1-B1./(Tb(i)+C1));
    P1 = P01.*x1(i);
    y1(i) = P1./P;
end
figure(5)
hold on
axis([0 1 min(Tb) max(Tb)])
plot(x1, Tb)
plot(y1, Tb)
xlabel('x1,y1')
ylabel('T [K]')

hold off

figure(6)
hold on
plot(x1,y1)
plot(x1,x1,'red')
xlabel('x')
ylabel('y')
% Molar flows 
F = sum(F_mol_cooler); xF=(sum((F_mol_cooler(1:3))))/F; xD=0.99; xB=0.01; %KÃ¤nt inflÃ¶de F o xF, Ã¶nskade sammansÃ¤ttningar xD(toppen) och xB=azeotrop(botten)
A=[1 1 F; xD xB xF*F]; Flows=rref(A); %rÃ¤knar ut toppflÃ¶de och bottenflÃ¶de med total samt komponentbalans
D=Flows(1,3); B=Flows(2,3);
q = 1;

%Calculate yF
Tbf = fsolve(@(T)find_Tbnonideal(T,xF,1-xF,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tbf+C1));
P1 = P01.*xF;
yF= P1./P; 

%Calculate R
LV = (xD - yF)/(xD - xF);
Rmin = LV/(1 - LV);
R = 1.5*Rmin;

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V;
Lbar = L + F;

%Initial temperature estimation and starting composition bottom of tower
Tstart = (tb1 + tb2)/2;
x(1) = xB;
Tb = fsolve(@(T)find_Tbnonideal(T,xB,1-xB,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tb+C1));
P1 = P01.*xB;
y(1)= P1./P;
% Botten
i = 1; 
while x<xF
    i = i + 1;
    x(i)=Vbar/Lbar*y(i-1) + B/Lbar*xB;
    y(i)=idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i));
    
end
% toppen
while y<xD
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F*xF);
    y(i)=idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i));
    i = i + 1;
end
n=i-1;
real=n/0.7;
% Torndimensioner

ts=0.45;
H=(real+1)*ts;

%Medelmolmassor
ML=x(1)*56.11+(1-x(1))*18.015; %  58.12  56.11
MV=y(1)*56.11+(1-y(1))*18.015;

%Temperatur i botten
Tbb = fsolve(@(T)find_Tbnonideal(T,x(1),1-x(1),1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
%Temperatur i toppen
Ttop = fsolve(@(T)find_Tbnonideal(T,xD,1-xD,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
%Medeldensiteter

MVrho=((ML/1000)*(P*133.322368))/(8.314*Tbb); %Kg/m^3
MLrho=(x(1))*588+(1-x(1))*997; %kg/mÂ³

FLV=(Lbar*ML)/(Vbar*MV)*sqrt(MVrho/MLrho);
CF=0.29; %flooding constant from diagram
sigma=72.8; %taken from some page water at 20Celcius
FST=(sigma/20)^0.2;
C=CF*FST;
Uf=C*sqrt(((MLrho-MVrho)/MVrho));
ada=0.1+(FLV-0.1)/9;
DT=sqrt((4*V*(MV/1000))/(0.8*(Uf/3.28)*pi*(1-ada)*MVrho));

% Totalkondensor och Ã¥terkokare
Hvap1=22.5*1000; %J/mol isobutene
Hvap2=44200; %J/mol vatten vid 20 grader C, 42911 vid 50 grader
Havgtop=xD*Hvap1+(1-xD)*Hvap2;
Havgbot=x(1)*Hvap1+(1-x(1))*Hvap2;

%condenser

Qc=D*(R+1)*Havgtop; %Joule/s

%reboiler
Qr=Vbar*Havgbot; %Joule/s
% COST DESTILLATIOn
% 2. COST: Distillation (Sieve trays=_s / Valve trays=_v / Bubble cap trays=_b)     %[m] 
%Constants
a_s=130;        a_v=210;        a_b=340;
b_s=440;        b_v=400;        b_b=640;
n_s=1.8;        n_v=1.9;        n_b=1.9;
%Calculations
%Pressure vesel
Pi=3.14;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)

Dim=DT;     %[m] diameter of container (1-3.5)
L=H;        %[m] length of container 
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=Pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=Pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;

%Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;     Cost_Year_B_h_k=a_h_k+b_h_k*Shell_mass.^n_h_k;     
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass_dist.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;     %[$] Cost calculations for different distillation trays

%Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_k=Cost_Year_B_h_k*(CEPCI_Year_A/CEPCI_Year_B);        
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

%Cost_Reactor_v_k=Cost_Year_A_v_k*10        %[SEK]
%Cost_Reactor_h_k=Cost_Year_A_h_k*10        %[SEK]
Cost_Dist_v_304_single=Cost_Year_A_v_304*10;       %[SEK]
Cost_Dist_h_304_single=Cost_Year_A_h_304*10;        %[SEK]

%Cost_Dist_v_304_five=Cost_Dist_v_304_single*ceil(real)
%Cost_Dist_h_304_five=Cost_Dist_h_304_single*ceil(real);

%Trays
S_dist=1.103;      %[m] diameter of trays
Cost_Year_B_s=a_s+b_s*S_dist.^n_s;     Cost_Year_B_v=a_v+b_v*S_dist.^n_v;     Cost_Year_B_b=a_b+b_b*S_dist.^n_b;     %[$] Cost calculations for different distillation trays

Cost_Year_A_s=Cost_Year_B_s*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_b=Cost_Year_B_b*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

Cost_Distillation_s=Cost_Year_A_s*10*ceil(real);        %[SEK]
Cost_Distillation_v=Cost_Year_A_v*10*ceil(real) ;       %[SEK]
Cost_Distillation_b=Cost_Year_A_b*10*ceil(real) ;       %[SEK]

Cost_Dist=Cost_Dist_v_304_single+Cost_Distillation_s;
fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*B*(xB),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*B*(xB),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*B*(xB),B*(1-xB),Tbb)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD),D*(1-xD),Ttop)
fprintf('Diameter of tower is %1.2f m and height is %1.2f m with %1.0f real trays. Cost of tower is %f SEK\n',Dim,H,ceil(real),Cost_Dist)
fprintf('Energy requirement in reboiler %1.2f MW. Energy requirement in condenser %1.2f MW\n\n\n',Qr*10^-6, Qc*10^-6)
%% Compressor 
fprintf('COMPRESSOR\n')
fprintf('Inflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD),D*(1-xD),Ttop)
format shortG

T_in= fsolve(@(T)find_Tbnonideal(T,xD,1-xD,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);


F_mol_compressor=[F_mol_cooler(1)/sum((F_mol_cooler(1:3)))*Flows(1,3)*xD;
       F_mol_cooler(2)/sum((F_mol_cooler(1:3)))*Flows(1,3)*xD;
       F_mol_cooler(3)/sum((F_mol_cooler(1:3)))*Flows(1,3)*xD;
       Flows(1,3)*(1-xD)];
   
F_mass_compressor=[F_mass(F_mol_compressor(1,1),1);
        F_mass(F_mol_compressor(2,1),2);
        F_mass(F_mol_compressor(3,1),3);
        F_mass(F_mol_compressor(4,1),4)];

CP_matrix=[Cp_new(T_in,1,1,1,1,1);
           Cp_new(T_in,2,2,2,2,2);
           Cp_new(T_in,3,3,3,3,3);
           Cp_new(T_in,4,4,4,4,4)];     %[J/Kg*K] Cp each component
  
% Calculations
C_matrix=F_mass_compressor.*CP_matrix;      
C_tot=sum(C_matrix);        %[W/K] Summan av m*cp för alla komponenter i flödet

P_in=1*10^5;     %[Pa] Ingående tryck till kompressorerna

P_ut=50*10^5;     %[Pa] Utgående tryck
eta_is=0.8;     %[] Isentropverkningsgrad
R=8.314;        %[J/mol*K]

M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

R_matrix=[R/M(1,1);
          R/M(2,1);
          R/M(3,1);
          R/M(4,1)];        %Butane-Butene-H2-H2O in J/kg*K
  
Cv_matrix=[CP_matrix(1,1)-R_matrix(1,1);
           CP_matrix(2,1)-R_matrix(2,1);
           CP_matrix(3,1)-R_matrix(3,1);
           CP_matrix(4,1)-R_matrix(4,1)];       %%Butane-Butene-H2-H2O
       
kappa_matrix=[CP_matrix(1,1)/Cv_matrix(1,1);
              CP_matrix(2,1)/Cv_matrix(2,1);
              CP_matrix(3,1)/Cv_matrix(3,1);
              CP_matrix(4,1)/Cv_matrix(4,1)];
          

kappa=sum(kappa_matrix)./4;     %!!!Dividedby number of components
       
% function [Wtot,Qkyltot,Akyltot,Tut]=kompressor(C_tot,kappa,Pin,Tin,Put,eta_is)
%Pressure increase per step
P_step = (P_ut/P_in)^(1/3);  %[]
%Temperature out for each step in isentrope compression
T_ut_is = T_in*P_step^((kappa-1)/kappa);  %[K] 
%Actual temperature out from every compression step.
T_ut = T_in + (T_ut_is-T_in)/eta_is; %[K] 
%Required compressor effekt for one compresson step.
W = C_tot*(T_ut-T_in); %[W] 
%Total required compression effect (3 steg).
W_tot = 3*W; %[W] 
%Required cooling effekt for 1 cooler inbetween steps 
Q_kyl = C_tot*(T_ut-T_in);%[W] 
%Total required cooling effect for all coolers between steps (2)
Q_kyl_tot = 2*Q_kyl; %[W] 
%Coolant temperature
T_kv = 14+273.15; %[K] 
%Maximal temperature that the coolant may be heated to
T_kv_max = 20+273.15; %[K] 
%Logarithmic average temperature difference
deltaTlm = ((T_in-T_kv)-(T_ut-T_kv_max))/log((T_in-T_kv)/(T_ut-T_kv_max)); %[]
%U-value for coolerbetween steps (gas-liquid)
Ukyl = 200; %[W/(m2K)] 
%Heat exchange unit for 1 cooler between steps
A_kyl = Q_kyl/(Ukyl*deltaTlm); %[m2] 
%Total heat exchange area for coolers 
A_kyl_tot = 2*A_kyl; %[m2] 

% Utdata
W_tot;        %[W] Totalt effektbehov för kompressionen
Q_kyl;      %[W] Kylbehov i mellankylare
A_kyl_tot;     %[m2] Total värmeväxlararea för mellankylare
T_ut;    %[K] Utgående temperatur
% 3. COST: Compressor      (Centrifugal compressor)        %[kW]
a=580000;
b=20000;
n=0.6;
S_comp=W_tot/1000;       %[kW]
Cost_Year_B=a+b*S_comp.^n;

Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_Compressor=Cost_Year_A*10;
fprintf('Outflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K and pressure %1.0f MPa\n', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD),D*(1-xD),T_ut,P_ut*10^-6)
fprintf('Cost of compressor is %f SEK\n\n\n',Cost_Compressor)

%% Cooler 2
format shortG
fprintf('COOLER 2\n')
fprintf('Inflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD),D*(1-xD),T_ut)
%Data
%Area calculations
C=0.8;      %Temperaturverkningsgrad from kurs PM
U=100;       %[W/m2] KVärmegenomgångstal from kurs PM

%Set up matrix for area calculations
Area_tot_1=zeros(10,1);

% Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler=F_mol_compressor;       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass_cooler=[F_mass(F_mol_cooler(1,1),1);
        F_mass(F_mol_cooler(2,1),2);
        F_mass(F_mol_cooler(3,1),3);
        F_mass(F_mol_cooler(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O

% Parameters for adjustment of loop
TH_in=T_ut;      %[K] Initial temperature of mixture
TH_ut_final=253;
TC_in=273-20;       %Inlet temperature of coolant
A_max=200;      %[m2] Maximum area per heat exchange unit
A_num_max=1000;       %Max number of heat exchange units

%Set up matrices for area calculations
m_1=zeros(A_num_max,1);       %Set up matrix for TC_out 
p_1=zeros(A_num_max,1);       %Set up matrix for TH_out
A_per_unit_1=zeros(A_num_max,1);       %Set up matrix for area per heat exchange unit
r_1=zeros(A_num_max,1);       %Set up matrix for number of heat exchange units

j=TH_in;      %Adjustment factor for TC_in
k=0;        %Adjustment factor for number of heat exchange units

% Loop calculating number of heat exchange units and area of units
while TH_in>TH_ut_final && k<A_num_max        %While loop that counts number or heat exchange units
k=k+1;
    while j>TH_ut_final       %While loop that iterates TC_in
    j=j-0.5;        %adjust for finer accuracy 
    
    i=0;
        while i<A_max       %While loop that tests area of each heat exchange unit
        i=i+0.01;        %adjust for finer accuracy    
        
        %Cp calculations mixture
        TH_ut=j;        %Outlet temperature
        TH_medel=(TH_in+TH_ut)/2;       

        CP_matrix=[Cp_new(TH_medel,1,1,1,1,1);
                   Cp_new(TH_medel,2,2,2,2,2);
                   Cp_new(TH_medel,3,3,3,3,3);
                   Cp_new(TH_medel,4,4,4,4,4)];     %[J/Kg*K] Cp each component
     
        %Temperature calculation coolant  
        TC_ut=TC_in+C*(TH_in-TH_ut);        %Outlet temperature of coolant     

        %Energy transfer calculations
        q_matrix=F_mass_cooler(:,1).*CP_matrix(:,1).*(TH_in-TH_ut);      %Energy transfer matrix
        q=sum(q_matrix);        %[J/s] Energy transfer

        %C_max/C_min calculations
        C_max_matrix=CP_matrix.*F_mass_cooler;   
        C_max=sum(C_max_matrix);
        C_min=C*C_max;
        
        %Epsilon calculation
        Epsilon=q/(C_min*(TH_in-TC_in));        %0 < Epsilon < 1
        
        %Function for checking area
        Fun=Epsilon-((1-exp(-(U.*i/C_min)*(1-(C_min/C_max))))/(1-(C_min/C_max)*exp(-(U.*i/C_min)*(1-C_min/C_max))));
        
            if Fun<0        %Checks to see if area is correct
                break;
            end
        end
            if i>A_max      %Checks if max area is exceeded 
                break;
            end
    end
A=i;
T=j;

%Outlet temperature for coolant 
TC_ut=TC_in+(TH_in-TH_ut)./C;

TH_in=T;
j=T;
m_1(k)=TC_ut;
p_1(k)=TH_in;
A_per_unit_1(k)=A;
r_1(k)=k;
end

% Output of heat exchange loop
Matrix_results=[r_1,m_1,p_1,A_per_unit_1];        %Matrix with results
column_names={'Number','TC_out','TH_out','Area'};        %Label for matrix with results 
Matrix_results( all(~Matrix_results,2),:)=[];
Heat_exchange_1=array2table(Matrix_results,'VariableNames',column_names)      %Matrix with results (labeled)

Area_tot_1=sum(A_per_unit_1);        %Total area for heat exchange
Area_tot_1=sprintf('%.f m2',Area_tot_1);        %Total area of all heat exchange units

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant
Cost_Year_B=zeros(10,1);

Shell_mass=A_per_unit(1,1);      %[m2] Area per unit
S_final=A_per_unit(k,1);        %[m2] Area of final unit

Cost_Year_B=(k-1)*(a_s+b_s*Shell_mass.^n_s)+(a_s+b_s*S_final.^n_s);      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_Cooler=Cost_Year_A*10;      %Conversion from Dollars to SEK'
fprintf('Total area needed is %s and total cost is %f SEK \n', Area_tot_1, Cost_Cooler)
fprintf('Outflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n\n\n', F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD),F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD),D*(1-xD),TH_ut_final)

%% Flash
%Antoine constants
Butanin=F_mol_cooler(1)/(sum((F_mol_cooler(1:3))))*D*(xD); Butenin=F_mol_cooler(2)/(sum((F_mol_cooler(1:3))))*D*(xD); H2in=F_mol_cooler(3)/(sum((F_mol_cooler(1:3))))*D*(xD); Win=D*(1-xD);
fprintf('FLASH\n')
fprintf('Inflow:%1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', Butanin,Butenin,H2in,Win,TH_ut_final)
A1 = [13.6333, 164.90, 3.19]; %hydorgen
A2 = [15.5381, 2032.73, -33.15]; %butan

%total pressure
P = 50*760;  %mmHg

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
figure(7)
plot(x1,T,'blue')
hold on
plot(y11,T,'red')
axis([0 1 min(T) max(T)])
xlabel('Fraction of "Hydrogen"')
ylabel('T [K]')
title('Temperature vs composition for hydrogen/"butane" system')
grid on
legend('Liquid phase composition', 'Vapor phase composition')

T=253; %temperature at which the flash will operate
fun=@(x1)(find_Tbflash(P,T,A1,A2,x1)-T);  %function to be solved
x11=fzero(fun,0.1); %solving the function to find the component fraction of the liquid phase
[T,y1]=find_Tbflash(P,T,A1,A2,x11); %using the found x-value to calculate the component fraction in the vapour phase

Fin =F_mol_cooler; %molar flows: Buthane|Buthene|Hydrogen|Water
Ftot=sum(Fin); 
z1in=Fin(3)/Ftot; %defining the inlet component fraction 
A=[1 1 Ftot;y1 x11 Ftot*z1in];%assembling the equations to be solved
flows=rref(A); %solving the equations by row reduction. First line gives flow for vapour, second gives liquid flow
hold on
line([z1in z1in],[0 T]); line([x11 y1], [T T]); line([x11 x11], [0 T], 'Linestyle', '--'); line([y1 y1], [0 T], 'Linestyle', '--') %drawing the flash operation into the diagram
legend('Liquid phase composition', 'Vapor phase composition')

% Dimensions

MV=y1*2*1.00784+(1-y1)*58.12; %average Mr of vapor g/mol
ML=x11*2*1.00784+(1-x11)*58.12; %average Mr of liquid g/mol

rhoV=((MV/1000)*(P*133.322368))/(8.314*T); %density of vapor kg/m3
rhoL=x11*70.85+(1-x11)*563; %density of liquid kg/m3

ut=0.07*sqrt((rhoL-rhoV)/rhoV); %m/s

L=flows(2,3)*(ML/1000)/rhoL; %liquid flow (mol/s)*(kg/mol)/(kg/m3)=m3/s
V=flows(1,3)*(MV/1000)/rhoV; %vapor flow (mol/s)*(kg/mol)/(kg/m3)=m3/s

DT=sqrt((4*V)/(pi*0.15*ut));

HL=(L*10*60)/((pi/4)*DT^2);

H=HL+1.5*DT;

%cost
Pi=3.14;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)

Dim=DT;     %[m] diameter of container (1-3.5)
L=H;        %[m] length of container 
S_max=137.9*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=50*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=Pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=Pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;

%Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;     Cost_Year_B_h_k=a_h_k+b_h_k*Shell_mass.^n_h_k;     
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass_dist.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;     %[$] Cost calculations for different distillation trays

%Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_k=Cost_Year_B_h_k*(CEPCI_Year_A/CEPCI_Year_B);        
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

%Cost_Reactor_v_k=Cost_Year_A_v_k*10        %[SEK]
%Cost_Reactor_h_k=Cost_Year_A_h_k*10        %[SEK]
Cost_Flash_v_304_single=Cost_Year_A_v_304*10;       %[SEK]
Cost_Flash_h_304_single=Cost_Year_A_h_304*10;        %[SEK]
fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n',Butanin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),Butenin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),x11*flows(2,3),Win/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),253)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', Butanin/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),Butenin/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),y1*flows(1,3),Win/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),253)
fprintf('Dimensions: diameter is %f and height is %f. The cost of the flash is %f\n\n\n',Dim,H,min(Cost_Flash_v_304_single,Cost_Flash_h_304_single))

%% Strypventil
Put=3.3*760; 
Tut=T; %antar att temperaturen ej förändras i strypventilen då vi rör oss inom underkyldvätska området

%% Destillation 2 butan-buten
%1 = buten, 2 = butan
butanin=Butanin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3); butenin=Butenin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3); h2in=x11*flows(2,3); win=Win/(Butanin+Butenin+Win)*(1-x11)*flows(2,3);
fprintf('DESTILLATION 2\n')
fprintf('Inflow:%1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', butanin,butenin,h2in,win,T)
clear x
clear y
clear i
clear n
%Wilson parameters
W12 = 0.48584; 
W21 = 1.64637;

%Antoine constants for degC, mmHg, log10

A2=15.7564; B2=2132.42; C2=-33.15 ;%buten
A1=15.5381; B1=2032.73; C1=-33.15;%butan
%total pressure
P =Put;  %mmHg

tb1=B1/(A1-log(P))-C1;
tb2=B2/(A2-log(P))-C2;

x1 = linspace(0,1,1000);
Tstart=(tb1+tb2)/2;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);

    %activity coefficients at x1
    gamma2 = exp(-log(x1(i)+W12*x2)+x2.*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
    gamma1 = exp(-log(x2+W21*x1(i))-x1(i).*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
   
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tbnonideal(T,x1(i),x2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = exp(A1-B1./(Tb(i)+C1));
    P1 = gamma1.*P01.*x1(i);
    y1(i) = P1./P;
end
figure(9)
hold on
axis([0 1 min(Tb) max(Tb)])
plot(x1, Tb)
plot(y1, Tb)
xlabel('x1,y1')
ylabel('T [K]')

hold off

figure(10)
hold on
plot(x1,y1)
plot(x1,x1,'red')
xlabel('x')
ylabel('y')
% Flows
F = butanin+butenin+h2in+win; xF=butanin/(butanin+butenin+h2in+win); xD=0.55; xB=0.1; %Känt inflöde F o xF, önskade sammansättningar xD(toppen) och xB=azeotrop(botten)
A=[1 1 F; xD xB xF*F]; Flows=rref(A); %räknar ut toppflöde och bottenflöde med total samt komponentbalans
D=Flows(1,3); B=Flows(2,3);


%Calculate yF
gamma2 = exp(-log(xF+W12*(1-xF))+(1-xF).*((W12./(xF+W12*(1-xF)))-(W21./(W21*xF+(1-xF)))));
gamma1 = exp(-log((1-xF)+W21*xF)-xF.*((W12./(xF+W12*(1-xF)))-(W21./(W21*xF+(1-xF)))));
Tbf = fsolve(@(T)find_Tbnonideal(T,xF,1-xF,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tbf+C1));
P1 = gamma1.*P01.*xF;
yF= P1./P; 
lambda=21.5*1000;
q=(95.21*(Tbf-T)+(lambda))/lambda;

%Calculate R
LV = (xD - yF)/(xD - xF);
Rmin = LV/(1 - LV);
R = 2*Rmin;

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V+F*(q-1);
Lbar = L + F*q;

%Initial temperature estimation and starting composition bottom of tower
Tstart = (tb1 + tb2)/2;
x(1) = xB;
gamma2 = exp(-log(xB+W12*(1-xB))+(1-xB).*((W12./(xB+W12*(1-xB)))-(W21./(W21*xB+(1-xB)))));
gamma1 = exp(-log((1-xB)+W21*xB)-xB.*((W12./(xB+W12*(1-xB)))-(W21./(W21*xB+(1-xB)))));
Tb = fsolve(@(T)find_Tbnonideal(T,xB,1-xB,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tb+C1));
P1 = gamma1.*P01.*xB;
y(1)= P1./P;
% Botten
i = 1; 
while x<xF
    i = i + 1;
    x(i)=Vbar/Lbar*y(i-1) + B/Lbar*xB;
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21);
    
end
% toppen
while y<xD
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F*xF);
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21);
    i = i + 1;
end
n=i-1;
real=n/0.7;
% Torndimensioner

ts=0.45;
H=(real+1)*ts;

%Medelmolmassor
ML=x(1)*58.12+(1-x(1))*56.11;
MV=y(1)*58.12+(1-y(1))*56.11;

%Temperatur i toppen
gamma2 = exp(-log(xD+W12*(1-xD))+(1-xD).*((W12./(xD+W12*(1-xD)))-(W21./(W21*xD+(1-xD)))));
gamma1 = exp(-log((1-xD)+W21*xD)-xD.*((W12./(xD+W12*(1-xD)))-(W21./(W21*xD+(1-xD)))));
Ttop = fsolve(@(T)find_Tbnonideal(T,x(1),1-x(1),gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);

%Medeldensiteter

MVrho=((ML/1000)*(P*133.322368))/(8.314*Tbb); %Kg/m^3
MLrho=(x(1))*563+(1-x(1))*588; %Kg/m^3

FLV=(Lbar*ML)/(Vbar*MV)*sqrt(MVrho/MLrho);
CF=0.28; %flooding constant from diagram
sigma=15.8; %taken from some sketchy page isobutene at 20Celcius
FST=(sigma/20)^0.2;
C=CF*FST;
Uf=C*sqrt(((MLrho-MVrho)/MVrho));
ada=0.1+(FLV-0.1)/9;
DT=sqrt((4*V*(MV/1000))/(0.8*(Uf/3.28)*pi*(1-ada)*MVrho));

% Totalkondensor och återkokare
Hvap1=21.5*1000; %J/mol isobutan
Hvap2=22.5*1000; %J/mol isobutene
Havgtop=xD*Hvap1+(1-xD)*Hvap2;
Havgbot=x(1)*Hvap1+(1-x(1))*Hvap2;

%condenser

Qc=D*(R+1)*Havgtop; %Joule/s

%reboiler
Qr=Vbar*Havgbot; %Joule/s
% 2. COST: Distillation (Sieve trays=_s / Valve trays=_v / Bubble cap trays=_b)     %[m] 
%Constants
a_s=130;        a_v=210;        a_b=340;
b_s=440;        b_v=400;        b_b=640;
n_s=1.8;        n_v=1.9;        n_b=1.9;
%Calculations
%Pressure vesel
Pi=3.14;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)

Dim=DT;     %[m] diameter of container (1-3.5)
L=H;        %[m] length of container 
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=Pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=Pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;

%Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;     Cost_Year_B_h_k=a_h_k+b_h_k*Shell_mass.^n_h_k;     
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass_dist.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;     %[$] Cost calculations for different distillation trays

%Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_k=Cost_Year_B_h_k*(CEPCI_Year_A/CEPCI_Year_B);        
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

%Cost_Reactor_v_k=Cost_Year_A_v_k*10        %[SEK]
%Cost_Reactor_h_k=Cost_Year_A_h_k*10        %[SEK]
Cost_Dist_v_304_single=Cost_Year_A_v_304*10;       %[SEK]
Cost_Dist_h_304_single=Cost_Year_A_h_304*10;        %[SEK]

%Cost_Dist_v_304_five=Cost_Dist_v_304_single*ceil(real)
%Cost_Dist_h_304_five=Cost_Dist_h_304_single*ceil(real);

%Trays
S_dist=1.103;      %[m] diameter of trays
Cost_Year_B_s=a_s+b_s*S_dist.^n_s;     Cost_Year_B_v=a_v+b_v*S_dist.^n_v;     Cost_Year_B_b=a_b+b_b*S_dist.^n_b;     %[$] Cost calculations for different distillation trays

Cost_Year_A_s=Cost_Year_B_s*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_b=Cost_Year_B_b*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

Cost_Distillation_s=Cost_Year_A_s*10*ceil(real);        %[SEK]
Cost_Distillation_v=Cost_Year_A_v*10*ceil(real) ;       %[SEK]
Cost_Distillation_b=Cost_Year_A_b*10*ceil(real) ;       %[SEK]

Cost_Dist=Cost_Dist_v_304_single+Cost_Distillation_s;

fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.3f K\n',butanin/(butanin+h2in)*xB*B,butenin/(butenin+win)*(1-xB)*B,h2in/(butanin+h2in)*(xB)*B,win/(butenin+win)*(1-xB)*B, Tb)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.3f K\n', butanin/(butanin+h2in)*D*xD, butenin/(butenin+win)*(1-xD)*D,h2in/(butanin+h2in)*(xD)*D,win/(butenin+win)*(1-xD)*D, Ttop)
fprintf('Diameter of tower is %1.2f m and height is %1.2f m with %1.0f real trays. Cost of tower is %f SEK\n',Dim,H,ceil(real),Cost_Dist)
fprintf('Energy requirement in reboiler %1.2f MW. Energy requirement in condenser %1.2f MW\n\n\n',Qr*10^-6, Qc*10^-6)




