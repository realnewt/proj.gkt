%Sep of "binary" sys. 
%1 = Buthane and H2/Buthene follows along, 2 = H2O

clear all

%Antoine constants for ln, K, mmHg
A1 = [15.5381, 2032.73, -33.15];
A2 = [18.3036, 3816.44, -46.13];

%Total pressure
P = 760;  %mmHg

%Molflöden: Buthane|Buthene|Hydrogen|Water

Fin = [11.5002, 42.4998, 42.4998, 540];


x1 = linspace(0, 1);
Tstart = 320;  %Temperature at which to start the search, in K

for i = 1:length(x1)
    x2 = 1 - x1(i);
   
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tb(T, x1(i), A1, A2, P), Tstart, options);
    
    P01 = exp(A1(1) - A1(2)./(Tb(i) + A1(3)));
    P1 = P01.*x1(i);
    y1(i) = P1./P;
    
    P011(i) = exp(A1(1) - A1(2)./(Tb(i) + A1(3)));

end

x1f = 0.1;
alfa = 0.5;
x1_op = 0:0.01:0.9;
y1_op = (x1f - (1 - alfa)*x1_op)/alfa;

Ustart = [Tstart, 0.5];
U = fsolve(@(U) find_flash(U, x1f, alfa, A1, A2, P), Ustart)