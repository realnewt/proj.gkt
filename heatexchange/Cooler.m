clear,clc;
format shortG

%% Data
%Molar flow
F_mol=[11.5002;
       42.9998;
       42.4998;
       540.0000];       %[mol/s] Butane-Butene-H2-H2O
      
%Area calculations
C=0.8;      %Temperaturverkningsgrad from kurs PM
U=50;       %[W/m2] KVärmegenomgångstal from kurs PM

%Cost calculations
a=32000;        %Cost constant
b=70;       %Cost constant
n=1.2;      %Equipment constant

%% Fetching and converting molar flow to mass flow
%Molar flow
F_mol=[11.5002;
       42.9998;
       42.4998;
       540.0000];       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass=[F_mass(F_mol(1,1),1);
        F_mass(F_mol(2,1),2);
        F_mass(F_mol(3,1),3);
        F_mass(F_mol(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O

%% Parameters for adjustment of loop
Cost=zeros(10,1);
Area_tot=zeros(10,1);
%for A_please = 1:5
TH_in=936;      %[K] Initial temperature of mixture
<<<<<<< Updated upstream
TH_ut_final=323;      %[K] Temperature after heat final exchange
A_max=350;      %[m2] Maximum area per heat exchange unit
=======
<<<<<<< Updated upstream
TH_ut_final=327;      %[K] Temperature after heat final exchange
A_max=5;      %[m2] Maximum area per heat exchange unit
=======
TH_ut_final=273.15+25.85;      %[K] Temperature after heat final exchange

A_max=650;      %[m2] Maximum area per heat exchange unit
>>>>>>> Stashed changes
>>>>>>> Stashed changes
A_num_max=1000;       %Max number of heat exchange units

%Set up matrices for area calculations
m=zeros(A_num_max,1);       %Set up matrix for TC_out 
p=zeros(A_num_max,1);       %Set up matrix for TH_out
A_per_unit=zeros(A_num_max,1);       %Set up matrix for area per heat exchange unit
r=zeros(A_num_max,1);       %Set up matrix for number of heat exchange units

j=936;      %Adjustment factor for TC_in
k=0;        %Adjustment factor for number of heat exchange units

%% Loop calculating number of heat exchange units and area of units
while TH_in>TH_ut_final && k<A_num_max        %While loop that counts number or heat exchange units
k=k+1;
    while j>TH_ut_final       %While loop that iterates TC_in
    j=j-0.1;        %adjust for finer accuracy 
    
    i=0;
        while i<A_max       %While loop that tests area of each heat exchange unit
        i=i+0.1;        %adjust for finer accuracy    
        
        %Cp calculations mixture
        TH_ut=j;        %Outlet temperature
        TH_medel=(TH_in+TH_ut)/2;       

        CP_matrix=[Cp_new(TH_medel,1,1,1,1,1);
                   Cp_new(TH_medel,2,2,2,2,2);
                   Cp_new(TH_medel,3,3,3,3,3);
                   Cp_new(TH_medel,4,4,4,4,4)];     %[J/Kg*K] Cp each component
     
        %Temperature calculation coolant
        TC_in=273+14;       %Inlet temperature of coolant  
        TC_ut=TC_in+C*(TH_in-TH_ut);        %Outlet temperature of coolant     

        %Energy transfer calculations
        q_matrix=F_mass(:,1).*CP_matrix(:,1).*(TH_in-TH_ut);      %Energy transfer matrix
        q=sum(q_matrix);        %[J/s] Energy transfer

        %C_max/C_min calculations
        C_max_matrix=CP_matrix.*F_mass;   
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

%% Output of heat exchange loop
Matrix_results=[r,m,p,A_per_unit];        %Matrix with results
column_names={'Number','TC_out','TH_out','Area'};        %Label for matrix with results 
Matrix_results( all(~Matrix_results,2),:)=[];
Heat_exchange=array2table(Matrix_results,'VariableNames',column_names)      %Matrix with results (labeled)

<<<<<<< Updated upstream
Area_tot=sum(A_per_unit);        %Total area for heat exchange
Area_tot=sprintf('%.f m2',Area_tot)
=======
Area_tot=sum(A_per_unit);        %Total area for heat exchange !!!
%Area_tot=sprintf('%.f m2',Area_tot)
>>>>>>> Stashed changes

%% Cost of heat exchange units

S=A_per_unit(1,1);      %[m2] Area per unit
S_final=A_per_unit(k,1);        %[m2] Area of final unit

<<<<<<< Updated upstream
Cost=(k-1)*(a+b*S.^n)+(a+b*S_final.^n);      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
Cost=CommaFormat(Cost); 
Cost=sprintf('$%s',Cost)


=======
Cost=(k-1)*(a+b*S.^n)+(a+b*S_final.^n);     %!!! %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
%Cost=CommaFormat(Cost); 
%Cost=sprintf('$%s',Cost);

%plot(Area_tot,Cost)
%end

%Area_tot
%Cost
%r;
%x=linspace(1,5);
%hold on
%plot(Area_tot,x)
%plot(Cost,x)
%hold off
>>>>>>> Stashed changes


