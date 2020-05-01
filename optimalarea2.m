function W=optimalarea2(U,Ka,Cmin,Cmax,Thin,Tcin,tdrift,beta,A0)


f =@(A)((Cmin^2)*exp(-(U*A/Cmin)*((-Cmin+Cmax)/Cmax))-2*Cmin*Cmax*exp(-(U*A/Cmin)*((-Cmin+Cmax)/Cmax))+(Cmax^2)*exp(-(U*A/Cmin)*((-Cmin+Cmax)/Cmax)))/(Cmax-Cmin*exp(-(U*A/Cmin)*((-Cmin+Cmax)/Cmax)))^2;
H=Ka/(U*(Thin-Tcin)*tdrift*beta);

L=@(A)f(A)-H;
options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt')
W=fsolve(L,A0, options);

end

