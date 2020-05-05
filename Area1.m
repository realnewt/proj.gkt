function a=Area1(E,c,Cmin,U)

%ntu=@(A)U*A/Cmin;
f=@(A)(1-exp(-(U.*A/Cmin)*(1-c)))/(1-(c*exp(-(U.*A/Cmin)*(1-c))))
options = optimoptions('fsolve','Algorithm',' Levenberg-Marquardt ')

sef=@(A)f(A)-E
a=fsolve(sef,2000,options)
end


%    