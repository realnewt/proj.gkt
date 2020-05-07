clear,clc;
format shortG

%% Data
T_in=273.15+400;

F_mol=[11.5002;
       42.9998;
       42.4998];
       %540.0000];
   
F_mass=[F_mass(F_mol(1,1),1);
        F_mass(F_mol(2,1),2);
        F_mass(F_mol(3,1),3)];
        %F_mass(F_mol(4,1),4)];

CP_matrix=[Cp_new(T_in,1,1,1,1,1);
           Cp_new(T_in,2,2,2,2,2);
           Cp_new(T_in,3,3,3,3,3)];
           %Cp_new(T_in,4,4,4,4,4)]
  
%% Calculations
C_matrix=F_mass.*CP_matrix;      
C_tot=sum(C_matrix);        %[W/K] Summan av m*cp f�r alla komponenter i fl�det

Pin=1*10^5;     %[Pa] Ing�ende tryck till kompressorerna
Tin=273.15+50;      %[K] Ing�ende temperatur
Put=3*10^5;     %[Pa] Utg�ende tryck
eta_is=0.8;     %[] Isentropverkningsgrad
R=8.314;        %[J/mol*K]

M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3];
   %18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

R_matrix=[R/M(1,1);
          R/M(2,1);
          R/M(3,1)];
          %R/M(4,1)];        %Butane-Butene-H2-H2O
  
Cv_matrix=[CP_matrix(1,1)-R_matrix(1,1);
           CP_matrix(1,1)-R_matrix(2,1);
           CP_matrix(1,1)-R_matrix(3,1)];
           %CP_matrix(1,1)-R_matrix(4,1)];       %%Butane-Butene-H2-H2O
       
kappa_matrix=[CP_matrix(1,1)/Cv_matrix(1,1);
              CP_matrix(2,1)/Cv_matrix(2,1);
              CP_matrix(3,1)/Cv_matrix(3,1)];
              %CP_matrix(4,1)/Cv_matrix(4,1)];
          
kappa=sum(kappa_matrix)./3;     %!!!Dividedby number of columns
        
%% function [Wtot,Qkyltot,Akyltot,Tut]=kompressor(C_tot,kappa,Pin,Tin,Put,eta_is)
%Tryck�kning per steg.
P_step = (Put/Pin)^(1/3);  %[]
%Temperatur ut fr�n varje kompressorsteg f�r isentrop kompression.
Tut_is = Tin*P_step^((kappa-1)/kappa);  %[K] 
%Verklig temperatur ut fr�n varje kompressorsteg.
Tut = Tin + (Tut_is-Tin)/eta_is; %[K] 
%Erforderlig kompressoreffekt f�r ett kompressorsteg.
W = C_tot*(Tut-Tin); %[W] 
%Total erforderlig kompressoreffekt (3 steg).
Wtot = 3*W; %[W] 
%Erforderlig kyleffekt i 1 mellankylare
Qkyl = C_tot*(Tut-Tin);%[W] 
%Total erforderlig kyleffekt i mellankylare (2 st)
Qkyltot = 2*Qkyl; %[W] 
%Kylvattnets temperatur.
Tkv = 14+273.15; %[K] 
%Maximal temperatur som kylvattnet f�r v�rmas till
Tkvmax = 20+273.15; %[K] 
%Logaritmisk medeltemperaturdifferens.
deltaTlm = ((Tin-Tkv)-(Tut-Tkvmax))/log((Tin-Tkv)/(Tut-Tkvmax)); %[]
%U-v�rde f�r mellankylare (gas-v�tska)
Ukyl = 200; %[W/(m2K)] 
%V�rmev�xlararea f�r 1 mellankylare
Akyl = Qkyl/(Ukyl*deltaTlm); %[m2] 
%Total v�rmev�xlararea f�r mellankylarna.
Akyltot = 2*Akyl; %[m2] 
%end

%% Utdata
Wtot        %[W] Totalt effektbehov f�r kompressionen
Qkyl        %[W] Kylbehov i mellankylare
Akyltot     %[m2] Total v�rmev�xlararea f�r mellankylare
Tut     %[K] Utg�ende temperatur





