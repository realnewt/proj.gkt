function Cp=Cp_new(T,a,b,c,d,m)      %Calculates Cp for each component

%Molar mass
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

%Specific heat capacities 
Cp_initial=[1.39 0.3847 -1.846e-04 2.895e-08;
            16.05 0.2804 -1.091e-04 9.098e-09;
            27.14 0.009274 -1.3813e-05 7.645e-09;
            32.24 0.001924 1.055e-05 -3.596e-09];       %Butane-Butene-H2-H2O in J/mol/K

Cp=(Cp_initial(a,1)+Cp_initial(b,2).*T+Cp_initial(c,3).*T^2+Cp_initial(d,4).*T^3)./M(m,1);   %based on c=(a+b*T+c*T^2+d*T^3)/M;        
end

