clear,clc;
format shortG

%%Konstanter
% Molmassa f�r komponenterna
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];                                                         %Butan-Buten-H2-H2O i kg/mol

%Specifika v�rmekapaciteterna f�r komponenterna
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09;
    32.24 0.001924 1.055e-05 -3.596e-09];                                   %Butan-Buten-H2-H2O i J/mol/K

C=0.8;                                                                      %Temperaturverkningsgrad f�s fr�n kurs PM
U=50;                                                                       %V�rmegenomg�ngstal fr�n kurs PM f�r gas n�ra atmosf�rtryck i W/m2K
%%V�rmev�xling fl�de 12
%molfl�de f�r fl�de 12 (varmt)
F_mol=[11.5002;
       42.9998;
       42.4998;
       540.0000];                                                           %Butan-Buten-H2-H2O i mol/s

%Butan-Buten-H2-H2O i kg/s       
F_massa=(F_mol(:,1).*M(:,1));                                               %Butan-Buten-H2-H2O i kg/s 
    
%function A=Please(Ji,j,k,l,J,I)
%%Kod f�r att ber�kna area av v�rmev�xlare
%Temp f�r fl�de 12
TH_in=936;                  %Initial temp
J=327;                          %Justerar initil temp (kanske ger litet fel om inte lika med slut temp T)
I=200;                         %Justerar till�ten area
l=10;
j=936;
k=0;

m=zeros(l,1);
n=zeros(l,1);
o=zeros(l,1);

while TH_in>J && k<l                 %Antal v�rmev�xlare
k=k+1;
while j>J                 %Slut temp (change at start of section)
    j=j-0.1;                %justera finhet f�r noggranare svar
    
i=0;
while i<I               %Max area of v�rmev�xlare
    i=i+0.1;                %justera finhet f�r noggranare svar    

TH_ut=j;                    %Slut temp
TH_medel=(TH_in+TH_ut)/2;                                                   %Temp i Kelvin
%Cp f�r componenter i fl�de 12

CP_flow_12=[Cp(TH_medel,CP(1,1),CP(1,2),CP(1,3),CP(1,4),M(1));
            Cp(TH_medel,CP(2,1),CP(2,2),CP(2,3),CP(2,4),M(2));
            Cp(TH_medel,CP(3,1),CP(3,2),CP(3,3),CP(3,4),M(3));
            Cp(TH_medel,CP(4,1),CP(4,2),CP(4,3),CP(4,4),M(4))];             %CP for flow 12 i J/Kg*K

%Cp och temp f�r kalt fl�de (vatten)
TC_in=273-20;                                                              %TC_ut=1219.3; %fr�n temp f�r utfl�de kallt  
TC_ut=TC_in+C*(TH_in-TH_ut);                                                %Temp av utfl�det av kalt fl�de i K    
TC_medel=(TC_in+TC_ut)/2;                                                   %Medel temp av utfl�det av kalt fl�de i K  

%q f�r v�rmev�xlare efter fl�de 12
q_matrix=F_massa(:,1).*CP_flow_12(:,1).*(TH_in-TH_ut);                      %Energi f�r upphettning
q=sum(q_matrix);                                                            %J/s

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
if i>I                    %Area av varje v�rmev�lare (change at start of sektion)
    break;
end
end
A=i;
T=j;

%Temp f�r utfl�de kallt
TC_ut=TC_in+(TH_in-TH_ut)./C;

TH_in=T;
j=T;
m(k)=TC_ut;
n(k)=TH_in;
o(k)=A;
end
yay=[n,o]
sum(o)

%A_fsolve=Area1(Epsilon,C_min,C_max,U)





