clear,clc;
format shortG

%%Konstanter
% Molmassa för komponenterna
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];                                 %Butan-Buten-H2-H2O i kg/mol

%Specifika värmekapaciteterna för komponenterna
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09;
    32.24 0.001924 1.055e-05 -3.596e-09];           %Butan-Buten-H2-H2O i J/mol/K

C=0.8;                                              %Temperaturverkningsgrad fås från kurs PM
U=50;
%%Värmeväxling flöde 12
%molflöde för flöde 12 (varmt)
F_mol=[11.5002;
       42.9998;
       42.4998;
       540.0000];                                   %Butan-Buten-H2-H2O i mol/s

%Butan-Buten-H2-H2O i kg/s       
F_massa=[F_mol(:,1).*M(:,1)];                        %Butan-Buten-H2-H2O i kg/s 
     
%Temp för flöde 12 
j=935;
while j>327
    j=j-0.1;
    
i=0;
while i<10000
    i=i+0.1;
TH_in=936;  TH_ut=j; TH_medel=(TH_in+TH_ut)/2;    %Temp i Kelvin

%Cp för componenter i flöde 12

CP_flow_12=[Cp(TH_medel,CP(1,1),CP(1,2),CP(1,3),CP(1,4),M(1));
            Cp(TH_medel,CP(2,1),CP(2,2),CP(2,3),CP(2,4),M(2));
            Cp(TH_medel,CP(3,1),CP(3,2),CP(3,3),CP(3,4),M(3));
            Cp(TH_medel,CP(4,1),CP(4,2),CP(4,3),CP(4,4),M(4))];     %CP for flow 12 i J/Kg*K

%Cp och temp för kalt flöde (vatten)
TC_in=273+180;                                                      %TC_ut=1219.3; %från temp för utflöde kallt  
TC_ut=TC_in+C*(TH_in-TH_ut);                                        %Temp av utflödet av kalt flöde i K    
TC_medel=(TC_in+TC_ut)/2;                                           %Medel temp av utflödet av kalt flöde i K  

%q för värmeväxlare efter flöde 12
q_matrix=F_massa(:,1).*CP_flow_12(:,1).*(TH_in-TH_ut);               %Energi för upphettning
q=sum(q_matrix);                                                     %J/s

%Temp för utflöde kallt
%TC_ut_2=TC_in+(TH_in-TH_ut)/C
C_max_matrix=CP_flow_12.*F_massa;   
C_max=sum(C_max_matrix);
C_min=C*C_max;

Epsilon=q/(C_min*(TH_in-TC_in));

Fun=Epsilon-((1-exp(-(U.*i/C_min)*(1-(C_min/C_max))))/(1-(C_min/C_max)*exp(-(U.*i/C_min)*(1-C_min/C_max))));
if Fun<0
   break;
end
end

if i>1000
    break;
end
end
clc
A=i
T=j

A_fsolve=Area1(Epsilon,C_min,C_max,U)





