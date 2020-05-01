%Kompressorberäkningar för GKT-projekt
%Gäller för alkandehydrerings-projekten samt metanolprojektet.
%Tänk på att endast gaser tryckhöjs i kompressorer. För vätskor används
%pumpar.
%Skapad av: Elin Göransson, 2008-04-07
%
%Beräkningar bygger på antaganden om adiabatisk kompression och omräkning
%med isentropverkningsgrad för att få verkligt effektbehov.
%
%Kompressionen delas upp i tre steg med mellankylning pga den stora
%tryckökningen. Uppdelningen sker så att effektbehovet blir samma i varje
%steg, och kylningen emellan utformas så att man får samma temperatur in
%till varje kompressor.
%
%Funktionen:
%[Wtot,Qkyl,Akyltot,Tut]=kompressor(Ctot,kappa,Pin,Tin,Put,eta_is)
%
%Indata:
%Ctot [W/K] = Summan av m*cp för alla komponenter i flödet, där m är
%             flödet i kg/s (alt mol/s) och cp är medelvärmekapaciviteten
%             över temperaturintervallet i kompressorn i J/(kgK) (alt J/(molK)).
%kappa []   = Kappatalet (viktat medelvärde av kappa för de olika
%             komponenterna)
%Pin [Pa]   = Ingående tryck till kompressorerna
%Tin [K]    = Ingående temperatur
%Put [Pa]   = Utgående tryck
%eta_is []  = Isentropverkningsgrad
%
%Utdata:
%Wtot [W]       = Totalt effektbehov för kompressionen.
%Qkyl [W]       = Kylbehov i mellankylare.
%Akyltot [m2]   = Total värmeväxlararea för mellankylare.
%Tut [K]        = Utgående temperatur.

function [Wtot,Qkyltot,Akyltot,Tut]=kompressor(Ctot,kappa,Pin,Tin,Put,eta_is)

%Tryckökning per steg.
P_step = (Put/Pin)^(1/3);  %[]
%Temperatur ut från varje kompressorsteg för isentrop kompression.
Tut_is = Tin*P_step^((kappa-1)/kappa);  %[K] 
%Verklig temperatur ut från varje kompressorsteg.
Tut = Tin + (Tut_is-Tin)/eta_is; %[K] 
%Erforderlig kompressoreffekt för ett kompressorsteg.
W = Ctot*(Tut-Tin); %[W] 
%Total erforderlig kompressoreffekt (3 steg).
Wtot = 3*W; %[W] 
%Erforderlig kyleffekt i 1 mellankylare
Qkyl = Ctot*(Tut-Tin);%[W] 
%Total erforderlig kyleffekt i mellankylare (2 st)
Qkyltot = 2*Qkyl; %[W] 
%Kylvattnets temperatur.
Tkv = 14+273.15; %[K] 
%Maximal temperatur som kylvattnet får värmas till
Tkvmax = 20+273.15; %[K] 
%Logaritmisk medeltemperaturdifferens.
deltaTlm = ((Tin-Tkv)-(Tut-Tkvmax))/log((Tin-Tkv)/(Tut-Tkvmax)); %[]
%U-värde för mellankylare (gas-vätska)
Ukyl = 200; %[W/(m2K)] 
%Värmeväxlararea för 1 mellankylare
Akyl = Qkyl/(Ukyl*deltaTlm); %[m2] 
%Total värmeväxlararea för mellankylarna.
Akyltot = 2*Akyl; %[m2] 