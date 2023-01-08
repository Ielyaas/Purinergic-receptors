clear
close all
clc

%% PARAMETER FILE

% P2Y1 parameters
param.A1 = 0; param.K1 = 0.05; param.K3 = 20;
param.k12 = 0.5; param.k_12 = 0.05; param.k14 = 0.05; param.k_14 = 0.5;

% P2Y2/4 parameters
param.A2 = 0.5; param.K2 = 0.05; param.K4 = 20;
param.k22 = 0.5; param.k_22 = 0.05; param.k24 = 0.05; param.k_24 = 0.5;

param.gamma = 5.5; param.delta = 0;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 2000; param.Ktau = 0.2;
param.Vplc = 0.0045; param.Kplc = 0.2; param.Vpkc = 0.2;
param.Kpkc = 0.3;

param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;

param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3;

param.L_t = 0.8; param.kL = 0.5; param.k_L = 0.008; param.k_r = 0.167;
param.Q_t1 = 0.887; param.Q_t2 = 1.087; param.kmin = 0.1; param.k_i = 0.1;

param.k_5P = 0.0001; param.k_3K = 0.05; param.K_3K = 0.4;

param.sigma1 = 1; param.sigma2 = 0.1;
param.omega1 = 0.5; param.omega2 = 0.8;
 
param.tau_c = 2; param.in = 1;
param.sw = 1; 

param.tau_L = 1;

param.lambda_1 = 15;               % Ktau
param.lambda_2 = 0.36;              % tau_max
param.lambda_3 = 1;                 % IP3 degradation
param.lambda_4 = 0.1;              % PKC-PLC interaction
param.lambda_5 = 75;                 % Influx PLC

%% ODE SOLVER

solver = @(x,t)Pur_rec(x,t,param);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[~,Y] = ode15s(solver,[0 5020],[1 1 0.0805296 3.861011 0.70395 0 0],opts);
[T,Y] = ode113(solver,[0 500],Y(end,:),opts);

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
% l_color = [0 0 1];
% r_color = [0 0.6 0];
% set(fig,'defaultAxesColorOrder',[l_color; r_color])
% yyaxis left
% plot(T,c,'LineWidth',2)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% % title('Halving V_{PLC}')
% % xlim([0 270])
% % set(gca,'yticklabel',[])
% hold on
% yyaxis right
% % plot(T,p,'color',[0 0.6 0],'LineWidth',4)
% xlabel('time (s)')
% ylabel('[IP_3] \muM')
% % set(gca,'yticklabel',[])
% hold on
% % plot(T,Jin,'k','LineWidth',2)
% % hold on
% % plot(T,Jpm,'c','LineWidth',2)
% % h = vline2(250,'r--');
% ax=gca;
% set(ax,'Linewidth',4)
% ax.FontSize=20;
% box off
% % axis([xmin, xmax, ymin, ymax])
% hold off
% % set(gcf,'position',[10,10,1250,800]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'Testfig12','epsc')

figure(2)
plot(T,c,'k','LineWidth',3)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] (\muM)')
% title('Closed cell')
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=40;
box off
% axis([xmin, xmax, ymin, ymax])
hold off
set(gcf,'position',[10,10,1250,800]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Pur_slide4b','epsc')

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
Q_i = kplus_i*param.Q_t1./(kplus_i + param.kmin);
Q_c = kplus_c*param.Q_t2./(kplus_c + param.kmin);

Q_1 = param.sigma1.*Q_i + param.omega1.*Q_c;
Q_2 = param.sigma2.*Q_i + param.omega2.*Q_c;

% PLC activation
bal = ((R1*param.A1/(param.K1 + param.A1)) + ...
     (R2*param.A2/(param.K2 + param.A2))).*...
     (param.Vplc + param.lambda_4.*Q_i);
kr = bal.*(param.lambda_5*alpha0_til + c.^2./(c.^2 + param.Kplc^2));

Kc_tilde = param.Kc.*(1-param.sw.*Q_c);
tau_max_til = param.tau_max.*(1-param.lambda_2.*param.delta);
Ktau_til = param.Ktau.*(1 - param.lambda_1.*Q_i);
phi_c = c.^4./(c.^4 + Kc_tilde.^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = tau_max_til*Ktau_til.^4./(c.^4 + Ktau_til.^4);


beta = phi_c.*phi_p.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

% IP3 degradation
eta = param.k_3K/(param.k_3K + param.k_5P);
nu_deg = (param.lambda_3.*Q_i + 1).*(eta.*c.^2./(c.^2 + param.K_3K^2) - (1 - eta));
tau_p = 1/((param.k_3K + param.k_5P).*(1 + param.lambda_3.*Q_i));

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jpm = param.delta.*param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = (1-R1)*rho1_d - R1.*Q_1*rho1_p;
out(2) = (1-R2)*rho2_d - R2.*Q_2*rho2_p;
out(3) = (Jipr - Jserca + Jin - Jpm)./param.tau_c;
out(4) = Jin - Jpm;
out(5) = (h_inf - h)./tau;
out(6) = (L - nu_deg.*p)./tau_p;
% out(6) = (L - param.k_i.*p);
out(7) = (kr.*(param.L_t - L) - param.k_r.*L)./param.tau_L;
out = out';

end                                      