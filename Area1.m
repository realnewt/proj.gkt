function a=Area1(Epsilon,C_min,C_max,U)

f=@(A)((1-exp(-(U*A/C_min)*(1-(C_min/C_max))))/(1-(C_min/C_max)*exp(-(U*A/C_min)*(1-C_min/C_max))));
options=optimoptions('fsolve','Algorithm',' Levenberg-Marquardt ')

sef=@(A)f(A)-Epsilon;
a=fsolve(sef,0,options);
end


%    