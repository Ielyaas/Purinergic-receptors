clear
close all
clc

%% PARAMETER FILE

% P2Y1 parameters
param.A = 0.5; param.K1 = 0.05; param.K3 = 20;
param.k2 = 0.5; param.k_2 = 0.05; param.k4 = 0.05; param.k_4 = 0.5;

param.gamma = 5.5; param.delta = 2.5;
param.Kc = 0.2; param.Kp = 0.3; param.Kf = 40; param.Kb = 0.4;
param.Kh = 0.1; param.tau_max = 2000; param.Ktau = 0.1;
param.V_plc = 0.09; param.Kplc = 0.11; param.V_pkc = 0.1;
param.Kpkc = 0.3;

param.Vs = 0.9; param.Kbar = 1.5e-5; param.Ks = 0.2;

param.alpha0 = 0.003; param.alpha1 = 0.01; param.Kce = 14;
param.Vpm = 0.07; param.Kpm = 0.3;

param.L_t = 0.8; param.kL = 0.5; param.k_L = 0.008; param.k_r = 0.167;
param.Q_t1 = 0.887; param.Q_t2 = 1.087; param.kmin = 0.1; param.k_i = 2;

param.V_5p = 0.66; param.V_3k = 0.1; param.K_3k = 0.4;

param.sigma = 1; % param.sigma2 = 0.1;
param.omega = 1; % param.omega2 = 0.8;

param.tau_c = 2;
param.sw = 1; 
param.lambda1 = 0.05; param.lambda2 = 0.36;

%% ODE SOLVER

solver = @(x,t)Pur_rec(x,t,param);
opts = odeset('RelTol',1e-9,'AbsTol',1e-12);
[~,Y] = ode15s(solver,[0 500],[1 0.0805296 3.861011 0.70395 0],opts);
[T,Y] = ode113(solver,[0 200],Y(end,:),opts);

%% OUTPUT/Plots

R = Y(:,1);
c = Y(:,2);
ct = Y(:,3);
h = Y(:,4);
p = Y(:,5);


alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

% P2Y1
rho_d = param.k4/(1+param.K3) + param.k_2/(1+(1/param.K3));
rho_p = param.k_4*param.K1/(param.K1 + param.A)...
        + param.k2*param.A/(param.K1 + param.A);
    
ce = param.gamma.*(ct - c);

% PLC activation
kr = (R*param.A/(param.K1 + param.A)).*...
     param.V_plc.*(c.^2./(c.^2 + param.Kplc^2));

L = kr*param.L_t./(kr + param.k_r);
D = param.kL.*L/param.k_L;

% PKC activation
Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus1 = param.V_pkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4);
kplus2 = param.V_pkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q1 = kplus1*param.Q_t1./(kplus1 + param.kmin);
Q2 = kplus2*param.Q_t2./(kplus2 + param.kmin);

Ktau_til = param.Ktau.*(1-param.lambda1.*Q1);

nu_deg = (Q1 + 1).*(param.V_5p + ...
         param.V_3k.*c.^2./(c.^2 + param.K_3k^2));

fig = figure;
figure(1)
l_color = [0 0 1];
r_color = [0 0.6 0];
set(fig,'defaultAxesColorOrder',[l_color; r_color])
yyaxis left
plot(T,c,'LineWidth',4)
xlabel('time (s)')
ylabel('[Ca^{2+}_i] \muM')
hold on
yyaxis right
% plot(T,p,'color',[0 0.6 0],'LineWidth',4)
xlabel('time (s)')
ylabel('[IP_3] \muM')
hold on
% plot(T,Jin,'k','LineWidth',2)
% hold on
% plot(T,Jpm,'c','LineWidth',2)
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',4)
ax.FontSize=20;
box off
% axis([xmin, xmax, ymin, ymax])
hold off
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'wu_figure4d','epsc')

% figure(2)
% plot(T,ct,'LineWidth',3)

% figure(3)
% plot(T,R1,'LineWidth',3)
% hold on
% plot(T,R2,'LineWidth',3)

% figure(4)
% plot(h,c,'LineWidth',3)

% figure(5)
% plot(ct,c,'LineWidth',3)

% figure(6)
% plot(T,Q1,'LineWidth',3)
% hold on
% plot(T,Q2,'LineWidth',3)
% hold on
% plot(T,Jin,'k','LineWidth',2)
% legend('Q1','Q2','Jin')

[pks,locs] = findpeaks(c);
locs = locs(pks>0.3);
hold on
plot(T(locs),Y(locs,2),"*")

per_Y = Y(locs(end-1):locs(end),:);
per_T = T(locs(end-1):locs(end));

figure(2)
plot(per_T,per_Y)

for i=1:5
    per(:,i) = interp1(per_T,per_Y(:,i),...
    linspace(T(locs(end-1)),T(locs(end)),10000),'spline');
end

T_per = linspace(T(locs(end-1)),T(locs(end)),10000);

figure(3)
plot(T_per,per)

Z = [T_per',per];
fileID = fopen('p2y_per09.dat','w');
fprintf(fileID,'%8d %8d %8d %8d %8d %8d\n',Z');
fclose(fileID);

tilefigs
%% FUNCTION FILE

function out = Pur_rec(~,in,param)

R = in(1);
c = in(2);
ct = in(3);
h = in(4);
p = in(5);

% if t < 100
%     param.delta = 2.5;
% else
%     param.delta = 0;
% end

alpha0_til = param.delta*param.alpha0;
alpha1_til = param.delta*param.alpha1;

% P2Y1 
rho_d = param.k4/(1+param.K3) + param.k_2/(1+(1/param.K3));
rho_p = param.k_4*param.K1/(param.K1 + param.A)...
        + param.k2*param.A/(param.K1 + param.A);

    
ce = param.gamma.*(ct - c);

% PLC activation
kr = (R*param.A/(param.K1 + param.A)).* ...
     (param.V_plc.*c.^2./(c.^2 + param.Kplc^2));

L = kr*param.L_t./(kr + param.k_r);
D = param.kL.*L/param.k_L;

% PKC activation
Jin = alpha0_til + alpha1_til*param.Kce^4./(ce.^4 + param.Kce^4);
kplus1 = param.V_pkc.*D.*alpha1_til.*param.Kce^4./(ce.^4 + param.Kce^4);
kplus2 = param.V_pkc.*D.*c.^2./(c.^2 + param.Kpkc^2);
Q1 = kplus1*param.Q_t1./(kplus1 + param.kmin);
Q2 = kplus2*param.Q_t2./(kplus2 + param.kmin);

Q_tilde = param.sigma.*Q1 + param.omega.*Q2;

Kc_tilde = param.Kc.*(1-param.sw.*Q2);
Ktau_til = param.Ktau.*(1-param.lambda1.*Q1);
tau_max_til = param.tau_max.*(1-param.lambda2.*param.delta);

phi_c = c.^4./(c.^4 + Kc_tilde^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);
h_inf = param.Kh^4./(c.^4 + param.Kh^4);
tau = tau_max_til.*Ktau_til^4./(c.^4 + Ktau_til^4);

beta = phi_c.*phi_p.*h;
alpha = phi_p_down.*(1 - phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

nu_deg = (Q1 + 1).*(param.V_5p + ...
         param.V_3k.*c.^2./(c.^2 + param.K_3k^2));
     
Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);
Jpm = param.delta.*param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = (1-R)*rho_d - R.*Q_tilde*rho_p;
out(2) = (Jipr - Jserca + Jin - Jpm)./param.tau_c;
out(3) = Jin - Jpm;
out(4) = (h_inf - h)./tau;
out(5) = param.kL.*L - nu_deg.*p;
out = out';

end