clear;
At=pi*(0.23/1000)^2;
Ae=pi*(2.25/1000)^2;
P0=0.9e+5;
% TC=0:5:50;
% T0=273+TC;
Ti=273+20;
gamma=1.13;
R=188;
cp=1630;
k=0.017;
mu=0.00829/1000; %dyanmic viscosity at 30degC
g=9.81;
D=0.003175;

mdotc=At.*P0.*((gamma./(R.*300)).*(2./(gamma+1)).^((gamma+1)./(gamma-1))).^0.5;
mdot=(1:0.1:6).*mdotc;
Re=4.*mdot./(pi*D*mu);


Qdot=6;
T0=Ti+Qdot./(mdot.*cp)
%Qdot=mdot.*cp.*(T0-Ti);

Ma=4.338;
%mdot=6*At*P0*((gamma/(R*T0))^0.5)*((gamma+1)/2)^((gamma+1)*0.5/(1-gamma));
%mdotc=At*P0*((gamma/(R*T0))*(2/(gamma+1))^((gamma+1)/(gamma-1)))^0.5;

Pe=P0*(1+(gamma-1)/2*Ma^2)^(gamma/(1-gamma));
Te=T0.*(1+(gamma-1)/2*Ma^2)^(-1);
ve=Ma.*(gamma.*R.*Te).^0.5;
Th=mdotc.*ve+Pe*Ae;
Isp=Th./(mdotc.*g);


% Dmin=(mdot/(148.842))^0.5;
% D=0.003175;
% Dc=0.025;
% Pr=mu*cp/k;
% Re=4*mdot/(pi*D*mu);
% %Nu=0.023*(Re^0.85)*(Pr^0.4)*(D/Dc)^0.1;
% 
% f=0.316*(1/Re)^0.25;
% Nu=(f/8)*(Re-1000)*Pr/(1+12.7*(f/8)^0.5*(Pr^(2/3)-1));
% h=Nu*k/D;
% 
% q=Qdot*0.5;
% Tf=450;
% Tw=1.1*Tf;
% L=q/(pi*D*h*(Tw-Tf));

% %figure;
% plot(mdot.*1000000,Isp);
% grid on;
% title("Isp achieved vs Mass flow rate for T_i = 273 K");
% xlabel('Mass flow rate (mg/s)');
% ylabel('Isp (s)');  
% legend("P = 2 W","P = 4 W","P = 6 W","P = 8 W","P = 10 W");
% hold on;

% %figure;
% plot(mdot.*1000000,Isp);
% grid on;
% title("Isp achieved vs Mass flow rate for W = 6 W");
% xlabel('Mass flow rate (mg/s)');
% ylabel('Isp (s)');  
% legend("T_i=-30^oC","T_i=-20^oC","T_i=-10^oC","T_i=0^oC","T_i=10^oC","T_i=20^oC");
% hold on;


figure;
plot(mdot.*1000000,Re);
text(1.101990970211485e+04, 2.278065345137164e-04, '\downarrow Re_Critical')
grid on;
title("Reynolds no vs Mass flow rate");
xlabel('Mass flow rate (mg/s)');
ylabel('Reynolds Number');  
hold on;