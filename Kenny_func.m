
%% Cooler       %Run secton to get Heat Exchange table and Total Area
clear,clc;
format shortG

%Data
%Area calculations
C=0.8;      %Temperaturverkningsgrad from kurs PM
U=50;       %[W/m2] KVärmegenomgångstal from kurs PM

%Set up matrix for area calculations
Area_tot=zeros(10,1);

% Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler=[7.5549;
        27.92;
        27.92;
       10.171];       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass_cooler=[F_mass(F_mol_cooler(1,1),1);
        F_mass(F_mol_cooler(2,1),2);
        F_mass(F_mol_cooler(3,1),3);
        F_mass(F_mol_cooler(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O

% Parameters for adjustment of loop
TH_in=936;      %[K] Initial temperature of mixture
TH_ut_final=323;      %[K] Temperature after heat final exchange
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
    j=j-30;        %adjust for finer accuracy 
    
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
Heat_exchange=array2table(Matrix_results,'VariableNames',column_names);      %Matrix with results (labeled)

Area_tot=sum(A_per_unit);        %Total area for heat exchange
Area_tot=sprintf('%.f m2',Area_tot);        %Total area of all heat exchange units

%% Compressor       %Run section to get total required effect for compression
%clear,clc;
format shortG

T_in=290;

F_mol_compressor=[7.1328;
       26.36;
       42.4998;
       0.25591];
   
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
T_in=323;      %[K] Ingående temperatur
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
Q_kyl;       %[W] Kylbehov i mellankylare
A_kyl_tot;     %[m2] Total värmeväxlararea för mellankylare
T_ut;     %[K] Utgående temperatur

%% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]
a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

% Cost of heat exchange units (uses Cooler.m)
Cost_Year_B=zeros(10,1);

Shell_mass=A_per_unit(1,1);      %[m2] Area per unit
S_final=A_per_unit(k,1);        %[m2] Area of final unit

Cost_Year_B=(k-1)*(a_s+b_s*Shell_mass.^n_s)+(a_s+b_s*S_final.^n_s);      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_Cooler=Cost_Year_A*10;      %Conversion from Dollars to SEK

%% 2. COST: Distillation (Sieve trays=_s / Valve trays=_v / Bubble cap trays=_b)     %[m] 
%Constants
a_s=130;        a_v=210;        a_b=340;
b_s=440;        b_v=400;        b_b=640;
n_s=1.8;        n_v=1.9;        n_b=1.9;
%Calculations
S_dist=2.5;      %[m] diameter of trays
Cost_Year_B_s=a_s+b_s*S_dist.^n_s;     Cost_Year_B_v=a_v+b_v*S_dist.^n_v;     Cost_Year_B_b=a_b+b_b*S_dist.^n_b;     %[$] Cost calculations for different distillation trays

Cost_Year_A_s=Cost_Year_B_s*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_b=Cost_Year_B_b*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

Cost_Distillation_s=Cost_Year_A_s*10;        %[SEK]
Cost_Distillation_v=Cost_Year_A_v*10;        %[SEK]
Cost_Distillation_b=Cost_Year_A_b*10;        %[SEK]

%% 3. COST: Compressor      (Centrifugal compressor)        %[kW]
a=580000;
b=20000;
n=0.6;
S_comp=W_tot/1000;       %[kW]
Cost_Year_B=a+b*S_comp.^n;

Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_Compressor=Cost_Year_A*10;

%% 4. COST: Reactor (Vertikal, kolstål=_v_k / Horisontell, kolstål=_h_k / Vertikal, 304-ss=_v_304 / Horisontell, 304-ss=_h_304)     %[kg] 
%Constants
%a_v_k=11600;        a_h_k=210;      
%b_v_k=34;       b_h_k=400;      
%n_v_k=0.85;     n_h_k=1.9;     %Temperature is too high

Pi=3.14;
Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)


D=1;     %[m] diameter of container (1-3.5)
L=2*D;        %[m] length of container (2*diameter from kurs PM)
S_max=74.5*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*D)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=Pi*((D/2)^2)*L;      %[m^3]Volume of full tank
V_inner=Pi*(((D/2)-t)^2)*(L-t);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_304=17400;        a_h_304=12800;
b_v_304=79;        b_h_304=73;
n_v_304=0.85;        n_h_304=0.85;

%Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;     Cost_Year_B_h_k=a_h_k+b_h_k*Shell_mass.^n_h_k;     
Cost_Year_B_v_304=a_v_304+b_v_304*Shell_mass.^n_v_304;       Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass.^n_h_304;     %[$] Cost calculations for different distillation trays

%Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_k=Cost_Year_B_h_k*(CEPCI_Year_A/CEPCI_Year_B);        
Cost_Year_A_v_304=Cost_Year_B_v_304*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

%Cost_Reactor_v_k=Cost_Year_A_v_k*10        %[SEK]
%Cost_Reactor_h_k=Cost_Year_A_h_k*10        %[SEK]
Cost_Reactor_v_304_single=Cost_Year_A_v_304*10;        %[SEK]
Cost_Reactor_h_304_single=Cost_Year_A_h_304*10;        %[SEK]

Cost_Reactor_v_304_five=Cost_Reactor_v_304_single*5;
Cost_Reactor_h_304_five=Cost_Reactor_h_304_single*5;

Cost_Katalyst_v_304=Cost_Reactor_v_304_five*0.5;
Cost_Katalyst_h_304=Cost_Reactor_h_304_five*0.5;

%% 5. COST: Flash       (Same as distillation)       %[m]
%Constants: Uses same a, b, n, and CEPCI values as distillation

%Calculations
S_flash=2.5;      %[m] diameter of tray?

Cost_Year_B_s=a_s+b_s*S_flash.^n_s;     Cost_Year_B_v=a_v+b_v*S_flash.^n_v;     Cost_Year_B_b=a_b+b_b*S_flash.^n_b;     %[$] Cost calculations for different distillation trays
Cost_Year_A_s=Cost_Year_B_s*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        Cost_Year_A_b=Cost_Year_B_b*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation

Cost_Flash_s=Cost_Year_A_s*10;        %[SEK]
Cost_Flash_v=Cost_Year_A_v*10;        %[SEK]
Cost_Flash_b=Cost_Year_A_b*10;        %[SEK]

%% Total cost for plant
Langfactor=4;

%Storage for values
Cost_Cooler_1=1;        %Parameters set:
Cost_Cooler_2=1;        %Parameters set:

Cost_Distillation_1=1;      %Parameters set:
Cost_Distillation_2=1;      %Parameters set:

%Full Cost for 
Equipment_Cost=Cost_Cooler+Cost_Distillation_s+Cost_Compressor+Cost_Reactor_h_304_five+Cost_Katalyst_h_304;     %Missing flash, two coolers, two distillation 
Cost_set_up=Equipment_Cost*Langfactor


