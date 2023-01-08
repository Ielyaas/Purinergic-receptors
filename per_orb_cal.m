clear
close all
clc

%% PARAMETER FILE

% P2Y1 parameters
param.A1 = 0.5; param.K1 = 0.05; param.K3 = 20;
param.k12 = 0.5; param.k_12 = 0.05; param.k14 = 0.05; param.k_14 = 0.5;

% P2Y2/4 parameters
param.A2 = 0.5; param.K2 = 0.05; param.K4 = 20;
param.k22 = 0.5; param.k_22 = 0.05; param.k24 = 0.05; param.k_24 = 0.5;

param.gamma = 5.5; param.delta = 0;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 200; param.Ktau = 0.09;
param.Vplc = 0.04; param.Kplc = 0.11; param.Vpkc = 0.1;
param.Kpkc = 0.3;

param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;

param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3;

param.L_t = 0.8; param.kL = 0.5; param.k_L = 0.008; param.k_r = 0.167;
param.Q_t1 = 0.887; param.Q_t2 = 1.087; param.kmin = 0.1; param.k_i = 2;

param.sigma1 = 1; param.sigma2 = 0.1;
param.omega1 = 1; param.omega2 = 0.8;

param.tau_c = 2; param.in = 1;
param.sw = 1;
param.sig_plc = 10;

%% ODE SOLVER

solver = @(x,t)Pur_rec(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[~,Y] = ode15s(solver,[0 500],[1 1 0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(solver,[0 200],Y(end,:),opts);

%% OUTPUT/Plots

R1 = Y(:,1);
R2 = Y(:,2);
c = Y(:,3);
ct = Y(:,4);
h = Y(:,5);
p = Y(:,6);


alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

ce = param.gamma.*(ct - c);

% PLC activation
kr = ((R1*param.A1/(param.K1 + param.A1)) + ...
     (R2*param.A2/(param.K2 + param.A2))).*...
     param.Vplc.*(param.sig_plc.*alpha0_til +...
     (c.^2./(c.^2 + param.Kplc^2)));

L = kr*param.L_t./(kr + param.k_r);
D = param.kL.*L/param.k_L;

% PKC activation
Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus1 = param.in.*(param.Vpkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4));
kplus2 = param.Vpkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q1 = kplus1*param.Q_t1./(kplus1 + param.kmin);
Q2 = kplus2*param.Q_t2./(kplus2 + param.kmin);

% fig = figure;
% figure(1)
% l_color = [0 0 1];
% r_color = [0 0.6 0];
% set(fig,'defaultAxesColorOrder',[l_color; r_color])
% yyaxis left
% plot(T,c,'LineWidth',4)
% xlabel('time (s)')
% ylabel('[Ca^{2+}_i] \muM')
% hold on
% yyaxis right
% plot(T,p,'color',[0 0.6 0],'LineWidth',4)
% xlabel('time (s)')
% ylabel('[IP_3] \muM')
% hold on
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % axis([xmin, xmax, ymin, ymax])
% hold off
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'wu_figure4d','epsc')

figure(1)
plot(T,c,'LineWidth',3)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')

[pks,locs] = findpeaks(c);
locs = locs(pks>0.5);
hold on
plot(T(locs),Y(locs,3),"*")

per_Y = Y(locs(end-1):locs(end),:);
per_T = T(locs(end-1):locs(end));

figure(2)
plot(per_T,per_Y)

for i=1:6
    per(:,i) = interp1(per_T,per_Y(:,i),...
    linspace(T(locs(end-1)),T(locs(end)),10000),'spline');
end

T_per = linspace(T(locs(end-1)),T(locs(end)),10000);

T_per = T_per - T_per(1,1);

figure(3)
plot(T_per,per)

Z = [T_per',per];
fileID = fopen('bm_per04.dat','w');
fprintf(fileID,'%8d %8d %8d %8d %8d %8d %8d\n',Z');
fclose(fileID);


%% FUNCTION FILE

function out = Pur_rec(~,in,param)

R1 = in(1);
R2 = in(2);
c = in(3);
ct = in(4);
h = in(5);
p = in(6);

% if t < 100
%     param.delta = 2.5;
% else
%     param.delta = 0;
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

% PLC activation
kr = ((R1*param.A1/(param.K1 + param.A1)) + ...
     (R2*param.A2/(param.K2 + param.A2))).*...
     param.Vplc.*(param.sig_plc.*alpha0_til +...
     (c.^2./(c.^2 + param.Kplc^2)));

L = kr*param.L_t./(kr + param.k_r);
D = param.kL.*L/param.k_L;

% PKC activation
Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus1 = param.in.*(param.Vpkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4));
kplus2 = param.Vpkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q1 = kplus1*param.Q_t1./(kplus1 + param.kmin);
Q2 = kplus2*param.Q_t2./(kplus2 + param.kmin);

Q_tilde_1 = param.sigma1.*Q1 + param.omega1.*Q2;
Q_tilde_2 = param.sigma2.*Q1 + param.omega2.*Q2;

Kp_tilde = param.Kp.*(1-param.sw.*Q2);
 
phi_c = c.^4./(c.^4 + param.Kc^4);
phi_p = p.^2./(p.^2 + Kp_tilde.^2);
phi_p_down = Kp_tilde.^2./(p.^2 + Kp_tilde.^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = param.tau_max*param.Ktau^4./(c.^4 + param.Ktau^4);

beta = phi_c.*phi_p.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jpm = param.delta.*param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = (1-R1)*rho1_d - R1.*Q_tilde_1*rho1_p;
out(2) = (1-R2)*rho2_d - R2.*Q_tilde_2*rho2_p;
out(3) = (Jipr - Jserca + Jin - Jpm)./param.tau_c;
out(4) = Jin - Jpm;
out(5) = (h_inf - h)./tau;
out(6) = param.kL.*L - param.k_i.*p;
out = out';

end