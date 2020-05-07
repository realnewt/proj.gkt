%Binary dest. of buthene (1) and buthane (2)

clear, clc

P = 760; %mmHg
q = 0;
zfa = 0.8;   % totaltkonc av 1 i infl�det
xfa = 0.8; % konc. av 1 i v�tskeinfl�det
xd = 0.99;
xw = 0.1;
F = 250;                  %feed
R = 3;                    %�terfl�desf�rh�llande

%Wilson parametrar
W_Be_Ba=[0.48584 1.64637];


%Antoine constants
A_Be = [15.7564, 2132.42, -33.15]; %N�n av v�rdena avviker mot Leylas
A_Ba = [15.5381, 2032.73, -33.15];

% totalbalans + Komponent Ballans f�r att ber�kna str�marnas storlekar
W =  F*(zfa-xd)/(xw-xd);
D = F - W;
% balances over condenser.
L = R*D;
V = L+D;

l = L+q*F; % L-streck
v = V-(1-q)*F; % V-streck

%�terkokare
y0 = 5/3-5/(3+4.5*xw);    % Kolla upp senare
x11 = (v*y0+W*xw)/l;         % KB �terkokare

% Avdrivardel   Stegning till feedbotten:
x(1) = x11;
i = 1;

while x(i)<xfa
    i = i + 1;
    y(i) = 5/3-5/(3+4.5*x(i-1)); % Kolla i h�r
    x(i+1) = (v/l).*y(i)+(W/l)*xw;  % Komponent ballans f�r avdrivardelen
end

m = i+1;

% F�rst�rkare
while y(i)<xd
    x(i+1) = (V/L)*y(i)+(1/L)*(W*xw-F*zfa);   %Komponent ballans �ver f�rst�rkardelen
    i = i+1;
    y(i) = 5/3-5/(3+4.5*x(i)); %Jv samband
end

%v�rmebehov �terkokare
%Ber�kning av temperaturen i �terkokaren  (Bubbelpunkt)
Tstart = 100;
options = optimset('Display', 'off');
TK = fsolve(@(T)findTbForDest(T, xw, W_Be_Ba, A_Be, A_Ba, P), Tstart, options);

TK = TK + 273.15;

% Sammans�ttning ut ur �terkokare
x1 = y0;              % Komponent A
x2 = 1-y0;            % Komponent B

%Ber�kning av entalpier i gasfas
aAg = 0.69381e5;
bAg = 0.6752e1;
cAg = 0.13199;
aBg = 0.31596e5;
bBg = 0.15841e2;
cBg = 0.15429;

HA = aAg+bAg*TK+cAg*TK^2;
HB = aBg+bBg*TK+cBg*TK^2;
Hblandningg = HA*x1 + HB*x2;

%Entalpier v�tskefas
aAv = 0.19534e5;
bAv = 0.63711e2;
cAv = 0.12206;
aBv = -0.12588e5;
bBv = 0.14150e2;
cBv = 0.23130;

hA = aAv+bAv*TK+cAv*TK^2;
hB = aBv+bBv*TK+cBv*TK^2;
hblandningv = hA*x1 + hB*x2;

qw = v*(Hblandningg  - hblandningv)/3600;        % �ngbildningsv�rme i kW

disp(['�verf�rd effekt   ' num2str(round(qw)),'  kW'])
disp(['bottnar nedre del    ' num2str(m)])
disp(['bottnar �vre del     ' num2str(i-m)])
disp(['bottnar totalt       ' num2str(i)])

figure(1);
plot(1:i,y,'r');
hold on
plot(1:i,x);
legend('y1','x1');
ylabel('x1,y1');
xlabel('Botten nr (R�knat fr�n �terkokaren)');

% Minsta �terfl�desf�rh�llandet

yfa=5/3-5/(3+4.5*xfa);
Rmin=(xd-yfa)/(yfa-xfa);            % (LV)min = Rmin/(Rmin+1)
disp(['Minsta �terfl�desf�rh�lland:    ',num2str(Rmin)])