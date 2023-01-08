clear
close all
clc
%% PLOT INSTRUCTIONS

% Fig 2b Activation
% A2 = 0, Vplc = 0.033, lambda_5 = 1
% Fig 2d
% A1 = 0, Vplc = 0.015, lamda_5 = 65
% Fig 2e
% Vplc = 0.0055

% Fig 3/4 Downregulation
% Dri = 0.3, Drc = 0.8

% Fig 6 Inhibition
% Vpkc = 0

% Fig 9/10 Ca2+-free
% delta = 0, lambda_1 = 225, lambda_2 = 125 (Simulation is numerically unstable)
% Can also switch tau_max = 2000, Ktau = 0.2 instead of changing lambda_1/2

%% PARAMETER FILE

% P2Y1 parameters
param.A1 = 0.5; param.K1 = 0.05; param.K3 = 20;
param.k12 = 0.5; param.k_12 = 0.05; param.k14 = 0.05; param.k_14 = 0.5;

% P2Y2/4 parameters
param.A2 = 0; param.K2 = 0.05; param.K4 = 20;
param.k22 = 0.5; param.k_22 = 0.05; param.k24 = 0.05; param.k_24 = 0.5;

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.1;
param.Vplc = 0.03; param.Kplc = 0.11; param.Vpkc = 0.1;
param.Kpkc = 0.3;

param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;

param.alpha0 = 0.004; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3;

param.L_t = 0.8; param.kL = 0.5; param.k_L = 0.008; param.k_r = 0.167;
param.Q_t1 = 0.887; param.Q_t2 = 1.087; param.kmin = 0.1; param.k_i = 1;

param.sigma1 = 1; param.sigma2 = 0.1;
param.omega1 = 0.5; param.omega2 = 0.8;
 
param.tau_c = 2; param.in = 1; 

param.lambda_1 = 0;               % Ktau
param.lambda_2 = 0;              % tau_max
param.lambda_3 = 1;                 % IP3 degradation
param.lambda_4 = 0.15;              % PKC-PLC interaction
param.lambda_5 = 1;                 % Influx PLC
param.lambda_6 = 1;

% PKC-DR

param.DRi = 0.3;
param.DRc = 0.8;

% param.DRi = 1;
% param.DRc = 1;

param.Tscale = 0.1;

% PLC rate constants
param.Lrate_1 = 100;                % A1
param.Lrate_2 = 2;                  % A2

%% ODE SOLVER

solver = @(x,t)Pur_rec(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[~,Y] = ode15s(solver,[0 5250],[1 1 0.0805296 3.861011 0.70395 0 0],opts);
[T,Y] = ode113(solver,[0 3000],Y(end,:),opts);

%% OUTPUT/Plots

R1 = Y(:,1);
R2 = Y(:,2);
c = Y(:,3);
ct = Y(:,4);
h = Y(:,5);
p = Y(:,6);
L = Y(:,7);


% fig = figure;
% figure(1)
% l_color = [0 0.447 0.741];
% r_color = [0.466 0.674 0.188];
% set(fig,'defaultAxesColorOrder',[l_color; r_color])
% yyaxis left
% plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% hold on
% yyaxis right
% plot(T,p,'color',[0.466 0.674 0.188],'LineWidth',3)
% xlabel('time (s)')
% ylabel('[IP_3] \muM')
% % xline(2000,'r--','LineWidth',5)
% hold on
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=40;
% box off
% % axis([xmin, xmax, ymin, ymax])
% hold off
% set(gcf,'position',[10,10,1200,800]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'Fig3sa','epsc')

figure(2)
plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
xlabel('time (s)')
ylabel('[Ca^{2+}_{i}] (\muM)')
title('P2Y_1 + P2Y_{2/4}')
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=70;
box off
% axis([xmin, xmax, ymin, ymax])
hold off
set(gcf,'position',[10,10,1400,1000]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig3sb.png')

% tilefigs
%% FUNCTION FILE

function out = Pur_rec(~,in,param)

R1 = in(1);
R2 = in(2);
c = in(3);
ct = in(4);
h = in(5);
p = in(6);
L = in(7);

% if t < 2000
%     param.Vpkc = 0.1;
% else
%     param.Vpkc = 0;
% end

alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

% P2Y1
rho1_d = param.k14/(1+param.K3) + param.k_12/(1+(1/param.K3));
rho1_p = param.k_14*param.K1/(param.K1 + param.A1)...
        + param.k12*param.A1/(param.K1 + param.A1);

% P2Y2/4
rho2_d = param.k24/(1+param.K4) + param.k_22/(1+(1/param.K4));
rho2_p = param.k_24*param.K2/(param.K2 + param.A2)...
        + param.k22*param.A2/(param.K2 + param.A2);
    
ce = param.gamma.*(ct - c);

D = param.kL.*L/param.k_L;

% PKC activation
Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus_i = param.in.*(param.Vpkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4));
kplus_c = param.Vpkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q_i = param.DRi*kplus_i*param.Q_t1./(kplus_i + param.kmin);
Q_c = param.DRc*kplus_c*param.Q_t2./(kplus_c + param.kmin);

Q_1 = param.sigma1.*Q_i + param.omega1.*Q_c;
Q_2 = param.sigma2.*Q_i + param.omega2.*Q_c;

% PLC activation
bal = ((R1*param.A1/(param.K1 + param.A1)) + ...
     (R2*param.A2/(param.K2 + param.A2))).*...
     (param.Vplc + param.A2*param.lambda_4.*Q_i);
kr = bal.*(param.lambda_5*alpha0_til + c.^2./(c.^2 + param.Kplc^2));

Kc_tilde = param.Kc.*(1-param.lambda_6.*Q_c);
tau_max_til = param.tau_max.*(1-param.lambda_2.*alpha0_til);
Ktau_til = param.Ktau.*(1 - param.lambda_1.*alpha0_til);
phi_c = c.^4./(c.^4 + Kc_tilde.^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = tau_max_til*Ktau_til.^4./(c.^4 + Ktau_til.^4);


beta = phi_c.*phi_p.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));


Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jpm = param.delta.*param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = param.Tscale.*((1-R1)*rho1_d - R1.*Q_1*rho1_p);
out(2) = param.Tscale.*((1-R2)*rho2_d - R2.*Q_2*rho2_p);
out(3) = param.Tscale.*(Jipr - Jserca + Jin - Jpm)./param.tau_c;
out(4) = param.Tscale.*(Jin - Jpm);
out(5) = param.Tscale.*(h_inf - h)./tau;
out(6) = param.Tscale.*(L - param.k_i.*p);
out(7) = param.Tscale*(param.A1*param.Lrate_1 + param.A2*param.Lrate_2).*...
         (kr.*(param.L_t - L) - param.k_r.*L);
out = out';

end                                      