clear 
close all
clc

%% PARAMETER FILE

% P2Y1 parameters
param.A1 = 0.5; param.K1 = 0.05; param.K3 = 20;
param.k2 = 0.5; param.k_2 = 0.05; param.k4 = 0.05; param.k_4 = 0.5;

param.gamma = 5.5; 
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; 
param.Kplc = 0.11;
param.Kpkc = 0.3;

param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;

param.alpha0 = 0.0027; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3;

param.L_t = 0.8; param.kL = 0.5; param.k_L = 0.008; param.k_r = 0.167;
param.Q_t1 = 0.887; param.Q_t2 = 1.087; param.kmin = 0.1; param.k_i = 1;

% PKC phosphorylation parameters
param.sigma1 = 1;
param.omega1 = 0.5;

param.tau_c = 2; param.in = 1;
param.sw = 1;

% PKC effects
param.lambda_3 = 0;
param.lambda_5 = 1;

% Parameters for PKC inhibition (Change Vpkc)

param.lambda_1 = 0;               % Ktau
param.lambda_2 = 0;               % tau_max
param.Vplc = 0.033;
param.Vpkc = 0.1;
param.delta = 2.5;
param.tau_max = 200; param.Ktau = 0.1;

% Parameters for Ca2+-free (Change delta)

% param.lambda_1 = 0;               % Ktau
% param.lambda_2 = 0;                 % tau_max
% param.Vplc = 0.02;
% param.Vpkc = 0.1;
% param.delta = 0;
% param.tau_max = 2000; param.Ktau = 0.2;

% PKC-DR 

% param.DRi = 0.1;
% param.DRc = 0.8;

param.DRi = 1;
param.DRc = 1;

param.Tscale = 0.1;

%% ODE SOLVER

solver = @(x,t)P2Y1_test(x,t,param);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[~,Y] = ode15s(solver,[0 6000],[1 0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(solver,[0 5000],Y(end,:),opts);

%% OUTPUT/Plots

R1 = Y(:,1);
c = Y(:,2);
ct = Y(:,3);
h = Y(:,4);
p = Y(:,5);

alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

% P2Y1
rho_d = param.k4/(1+param.K3) + param.k_2/(1+(1/param.K3));
rho_p = param.k_4*param.K1/(param.K1 + param.A1)...
        + param.k2*param.A1/(param.K1 + param.A1);
    
ce = param.gamma.*(ct - c);

% PLC activation

kr = (R1*param.A1/(param.K1 + param.A1)).*param.Vplc.*...
     (param.lambda_5*alpha0_til + c.^2./(c.^2 + param.Kplc^2));
L = kr*param.L_t./(kr + param.k_r);

D = param.kL.*L/param.k_L;

% PKC activation

Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus_i = param.in.*(param.Vpkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4));
kplus_c = param.Vpkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q_i = param.DRi*kplus_i*param.Q_t1./(kplus_i + param.kmin);
Q_c = param.DRc*kplus_c*param.Q_t2./(kplus_c + param.kmin);

% 
fig = figure;
figure(1)
l_color = [0 0.447 0.741];
r_color = [0.466 0.674 0.188];
set(fig,'defaultAxesColorOrder',[l_color; r_color])
yyaxis left
plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
% title('Halving V_{PLC}')
% xlim([0 270])
% set(gca,'yticklabel',[])
hold on
yyaxis right
plot(T,p,'color',[0.466 0.674 0.188],'LineWidth',3)
xlabel('time (s)')
ylabel('[IP_3] \muM')
hold on
% plot(T,L,'r','LineWidth',3)
% set(gca,'yticklabel',[])
hold on
% plot(T,Jin,'k','LineWidth',2)
% hold on
% plot(T,Jpm,'c','LineWidth',2)
xline(2000,'r--','LineWidth',5)
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=40;
% set(gcf, 'color', 'none');
box off
% axis([xmin, xmax, ymin, ymax])
hold off
set(gcf,'position',[10,10,1200,800]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig6sb','epsc')

% figure(2)
% plot(T,c,'Color',[0 0.447 0.741],'LineWidth',3)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] (\muM)')
% % title('Closed cell')
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=70;
% box off
% % axis([xmin, xmax, ymin, ymax])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'Fig4sa','epsc')

% tilefigs

%% FUNCTION FILE

function out = P2Y1_test(t,in,param)

R1 = in(1);
c = in(2);
ct = in(3);
h = in(4);
p = in(5);

% if t < 2000
%     param.Vpkc = 0.1;
% else
%     param.Vpkc = 0;
% end
        
% if t < 500
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end
%         
alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

% P2Y1
rho_d = param.k4/(1+param.K3) + param.k_2/(1+(1/param.K3));
rho_p = param.k_4*param.K1/(param.K1 + param.A1)...
        + param.k2*param.A1/(param.K1 + param.A1);
    
ce = param.gamma.*(ct - c);

% PLC activation

kr = (R1*param.A1/(param.K1 + param.A1)).*param.Vplc.*...
     (param.lambda_5*alpha0_til + c.^2./(c.^2 + param.Kplc^2));
L = kr*param.L_t./(kr + param.k_r);

D = param.kL.*L/param.k_L;

% PKC activation

Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus_i = param.in.*(param.Vpkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4));
kplus_c = param.Vpkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q_i = param.DRi*kplus_i*param.Q_t1./(kplus_i + param.kmin);
Q_c = param.DRc*kplus_c*param.Q_t2./(kplus_c + param.kmin);

Q_1 = param.sigma1.*Q_i + param.omega1.*Q_c;

Kc_tilde = param.Kc.*(1-param.sw.*Q_c);
tau_max_til = param.tau_max.*(1-param.lambda_2.*alpha0_til);
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
% eta = param.k_3K/(param.k_3K + param.k_5P);
% nu_deg = (param.lambda_3.*Q_i + 1).*(eta.*c.^2./(c.^2 + param.K_3K^2) - (1 - eta));
% tau_p = 1/((param.k_3K + param.k_5P).*(1 + param.lambda_3.*Q_i));

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jpm = param.delta.*param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = param.Tscale.*((1-R1)*rho_d - R1.*Q_1*rho_p);
out(2) = param.Tscale.*(Jipr - Jserca + Jin - Jpm)./param.tau_c;
out(3) = param.Tscale.*(Jin - Jpm);
out(4) = param.Tscale.*(h_inf - h)./tau;
out(5) = param.Tscale.*(L - param.k_i.*p);
out = out';         

end 

