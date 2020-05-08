function a=Basic_data(F_mol,M,T_in,CP)      %Calculates Cp for each component

a=F_mol*M;

CP_matrix=[Cp(T_in,CP(1,1),CP(1,2),CP(1,3),CP(1,4),M(1));
           Cp(T_in,CP(2,1),CP(2,2),CP(2,3),CP(2,4),M(2));
           Cp(T_in,CP(3,1),CP(3,2),CP(3,3),CP(3,4),M(3));
           Cp(T_in,CP(4,1),CP(4,2),CP(4,3),CP(4,4),M(4))];
end
