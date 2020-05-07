%{
% A=butane B=butene H=hydrogen W=water
% matlabprogram to determine the number of reactors in the reaction 
%Butane --> Butene+H2 with assumed amount of katalyst and heating in
%between the reactors.
%}
clear, clc;
%clf
format short
e=0;
cat_tot=500; %kg cat
XA_start=0;
TCin=950; %constant start temp. for every reactor in kelvin
FA0=54; FB0=0.5; FW0=10*FA0; FH0=0; %molar flowrates in to the first reactor mol/s
C=0.8; U=50;
HR=116.3e3; %J/mol reaction enthalpy
P=1; %bar constant total pressure
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09;
    32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter J/mol/K
    
Tslut=zeros(1,5); %define a matrix just so matlab doesn't complain
XAut=zeros(1,5);
while XA_start-0.78<1e-4
%{
a while-loop to determine the amount of reactors needed 
to get the conversion sought after. the loop stops when converison hits a
specific level (0.95 in this case).
%}
[cat,Y]=ode15s(@(cat,Y) ode_func(cat, Y, HR, P, CP, FA0,FH0, FB0, FW0), [0 cat_tot], [XA_start TCin]); 
% differential equation solver where ode_func is the differential function
% at hand and the rest is the parameters and cat & Y matrixes. 

XA=Y(:,1); T=Y(:,2);


figure(1)
plot(cat,XA,'linewidth',1), hold on               % plot with conversion against catalyst mass.
xlabel('Mängd katalysator (kg)')
ylabel('XA')

figure(2)
set ( gca, 'xdir', 'reverse' )            %x axis get set to reverse number order (largest to smallest)
plot(T,XA,'linewidth',1), hold on                   %plot with conversion against temperature 
xlabel('Temperatur (K)')
ylabel('XA')



Tslut(e+1)=T(end); %end temperatures of each reactor for future use 
XA_start=max(XA);  %set the new start conversion for the new reactor as the end conversion of the reactor in this itteration
fa=(FA0)*(1-XA_start);
fb=FB0+FA0*XA_start; fh=FH0+FA0*XA_start; fw=FW0;
ftot=fa+fb+fh+fw;
XAut(e+1)=XA_start(end); 
e=e+1;                  %    counting the number of reactors
ff(e,:)=[fa fb fh fw];
leg(e,:)= "Reaktor "+ e;                                           %    matrix for the legend of the graphs
end
%{
figure(1)
legend(leg,'location','southeast') %legend with specific placement
figure(2)
legend(leg,'location','northeast') % -----------  I  I  ---------------
disp("reaktorer: "+e)                    % displaying the number of reactors.
%}
disp('grattis')

%% potentiell värmeväxlare
% Molmassa för komponenterna
M=[58.12*10^-3 ;
   56.1*10^-3 ;
   2*1.00784*10^-3 ;
   18.01528*10^-3];                             %Butan-Buten-H2-H2O i kg/mol
g=0;
for i=1:4
Fun=0.1;
disp('ny')
l=0;
while Fun>=0
    %  molflöde  + massa värme växlare 1
F_mol=[ff(1,i);
       ff(2,i);
       ff(3,i);
       ff(4,i)];
F_massa=[ff(1,i)*M(1);
             ff(2,i)*M(2);
             ff(3,i)*M(3);
             ff(4,i)*M(4);];                    %Butan-Buten-H2-H2O i kg/s  
     TCin=Tslut(i); TCut=950;   TCmedel=(TCin+TCut)/2;
cpBA=@(TCmedel)CP(1,1)+CP(1,2).*TCmedel+CP(1,3).*TCmedel.^2+CP(1,4).*TCmedel.^3; cpBA=cpBA(TCmedel);
cpBE=@(TCmedel)CP(2,1)+CP(2,2).*TCmedel+CP(2,3).*TCmedel.^2+CP(2,4).*TCmedel.^3; cpBE=cpBE(TCmedel);
cpH=@(TCmedel)CP(3,1)+CP(3,2).*TCmedel+CP(3,3).*TCmedel.^2+CP(3,4).*TCmedel.^3; cpH=cpH(TCmedel);
cpW=@(TCmedel)CP(4,1)+CP(4,2).*TCmedel+CP(4,3).*TCmedel.^2+CP(4,4).*TCmedel.^3; cpW=cpW(TCmedel);

cp=[cpBA
    cpBE
    cpH
    cpW];

QBA=F_massa(1)*cpBA*(TCut-TCin);
QBE=F_massa(2)*cpBE*(TCut-TCin);
QH=F_massa(3)*cpH*(TCut-TCin);
QW=F_massa(4)*cpW*(TCut-TCin);
Qtot=QBA+QBE+QH+QW;

%{
bestäm temperaturen för den värmande ångan???
huur?
%}
THin=1100;          
THut=THin+C*(THin-TCut);         
THmedel=(THin+THut)/2;

CC=sum(cp.*F_massa); Cmax=CC;
Cmin=C*Cmax;    %CH
%CH=cp(4)*F_massa(4)

l=l+1;

Epsilon=Qtot/(Cmin*(THin-TCin));
Fun=Epsilon-((1-exp(-(U.*l/Cmin)*(1-(Cmin/Cmax))))/(1-(Cmin/Cmax)*exp(-(U.*l/Cmin)*(1-Cmin/Cmax))));

end
l
g=g+l;
end % känns rätt såhär men lär vara fel.

%% ugn
clc
for i=1:4
  M=[58.12*10^-3 ;
   56.1*10^-3 ;
   2*1.00784*10^-3 ;
   18.01528*10^-3];

F_massa=[ff(1,i)*M(1);
             ff(2,i)*M(2);
             ff(3,i)*M(3);
             ff(4,i)*M(4);];
         
Tin=Tslut(i); Tut=950; Tmedel=(Tin+Tut)/2;
cpBA=@(Tmedel)CP(1,1)+CP(1,2).*Tmedel+CP(1,3).*Tmedel.^2+CP(1,4).*Tmedel.^3; cpBA=cpBA(Tmedel);
cpBE=@(Tmedel)CP(2,1)+CP(2,2).*Tmedel+CP(2,3).*Tmedel.^2+CP(2,4).*Tmedel.^3; cpBE=cpBE(Tmedel);
cpH=@(Tmedel)CP(3,1)+CP(3,2).*Tmedel+CP(3,3).*Tmedel.^2+CP(3,4).*Tmedel.^3; cpH=cpH(Tmedel);
cpW=@(Tmedel)CP(4,1)+CP(4,2).*Tmedel+CP(4,3).*Tmedel.^2+CP(4,4).*Tmedel.^3; cpW=cpW(Tmedel);

cp=[cpBA
    cpBE
    cpH
    cpW];

QBA=F_massa(1)*cpBA*(Tut-Tin);
QBE=F_massa(2)*cpBE*(Tut-Tin);
QH=F_massa(3)*cpH*(Tut-Tin);
QW=F_massa(4)*cpW*(Tut-Tin);
Qtot=QBA+QBE+QH+QW;
         
chi=0.8;
Qheat=chi*Qtot
end
