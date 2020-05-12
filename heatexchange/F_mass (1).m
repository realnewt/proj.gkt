function F_mass=F_mass(F_mol,m)

%Molar mass
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

F_mass=(F_mol.*M(m,1));       %Butane-Butene-H2-H2O in kg/s
end

