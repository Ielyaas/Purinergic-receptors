%% Experimental data for ATP

clear
close all
clc

%% DATA

[v,T,vT] = xlsread('ATP_data.xlsx');

%% PLOTS

% Ca response
t1 = v(:,1);
y1 = v(:,2);

% % ATP
figure(1)
plot(t1,y1,'k','LineWidth',4)
xlabel('time (s)')
ylabel({'Fura-2 ratio'})
hold on
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',6)
ax.FontSize=70;
box off
% yticks([0.9 1 1.1 1.2 1.3 1.4])
xlim([0 1600])
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
saveas(gcf,'Fig3a','epsc')