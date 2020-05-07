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
Tstart = 50 + 273.15;  %Temperature at which to start the search, in K

for i = 1:length(x1)
    x2 = 1 - x1(i);
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)findTbForFlash(T, x1(i), A1, A2, P), Tstart, options);
end

x1f = 0.1;
alfa = 0.5;

Ustart = [Tstart, 0.5];
U = fsolve(@(U) find_flash(U, x1f, alfa, A1, A2, P), Ustart)