%Sep of "binary" sys. 
%1 = H2/Buthane/Buthene, 2 = H2O

%Wilson parameters
W12 = ()/3;
W21 = ;

%Antoine constants for               degC, mmHg, log10
A1 = ;  B1 = ;  C1 = ;
A2 = ;  B2 = ;  C2 = ;

%total pressure
P = ;  %mmHg

x1 = linspace(0, 1);
Tstart=;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);

    %activity coefficients at x1
    gamma1 = exp(-log(x1(i)+W12*x2)+x2.*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
    gamma2 = exp(-log(x2+W21*x1(i))-x1(i).*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
   
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tb(T,x1(i),gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = 10.^(A1-B1./(Tb(i)+C1));
    P1 = gamma1.*P01.*x1(i);
    y1(i) = P1./P;
    
    P011(i) = 10.^(A1-B1./(Tb(i)+C1));

end

x1f = 0.5;
alfa = 0.5;
x1_op = 0:0.01:0.9;
y1_op = (x1f - (1-alfa)*x1_op)/alfa;

Ustart = [Tstart 0.5];
U = fsolve(@(U) find_flash(U, x1f, alfa, W12, W21, A1, A2, B1, B2, C1, C2, P), Ustart);