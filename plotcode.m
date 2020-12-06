clear;
At=pi*(0.23/1000)^2;
Ae=pi*(2.25/1000)^2;
P0=0.9e+5;
TC=0:5:50;
T0=273+TC;
Ti=273+20;
gamma=1.13;
R=188;
cp=1630;
k=0.017;
mu=0.00829/1000; %dyanmic viscosity at 30degC
g=9.81;

Ma=4.338;
%mdot=6*At*P0*((gamma/(R*T0))^0.5)*((gamma+1)/2)^((gamma+1)*0.5/(1-gamma));
%mdotc=At*P0*((gamma/(R*T0))*(2/(gamma+1))^((gamma+1)/(gamma-1)))^0.5;
mdotc=At.*P0.*((gamma./(R.*T0)).*(2./(gamma+1)).^((gamma+1)./(gamma-1))).^0.5;
mdot= 4*mdotc;
Pe=P0*(1+(gamma-1)/2*Ma^2)^(gamma/(1-gamma));
Te=T0.*(1+(gamma-1)/2*Ma^2)^(-1);
ve=Ma.*(gamma.*R.*Te).^0.5;
Th=mdotc.*ve+Pe*Ae;
Isp=Th./(mdotc.*g);

Qdot=mdot.*cp.*(T0-Ti);

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

%figure;
plot(Isp,Qdot);
title("Power required vs Isp")
xlabel('Specific Impulse (s)')
ylabel('Power Input (W)')
legend("T_i=-30^oC","T_i=-20^oC","T_i=-10^oC","T_i=0^oC","T_i=10^oC","T_i=20^oC")
hold on
