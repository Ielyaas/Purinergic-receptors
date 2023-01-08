%% Experimental data for figures

clear
close all
clc

%% DATA

[v,T,vT] = xlsread('Exp_data.xlsx');

%% PLOTS

% Ca response
t1 = v(:,1);
y1 = v(:,2);
t2 = v(:,5);
y2 = v(:,6);
t3 = v(:,9);
y3 = v(:,10);

% PKC inhibition
t4 = v(:,13);
y4 = v(:,14);
t5 = v(:,17);
y5 = v(:,18);

% PKC-DR
t6 = v(:,21);
y6 = v(:,22);
t7 = v(:,24);
y7 = v(:,25);

t8 = v(:,28);
y8 = v(:,29);
t9 = v(:,31);
y9 = v(:,32);

% Ca-free
t10 = v(:,35);
y10 = v(:,36);
t11 = v(:,38);
y11 = v(:,39);

t12 = v(:,42);
y12 = v(:,43);
t13 = v(:,45);
y13 = v(:,46);

% % ATP
% figure(1)
% plot(t1,y1,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% title('ATP')
% hold on
% % h = vline2(250,'r--');
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % yticks([0.9 1 1.1 1.2 1.3 1.4])
% xlim([0 1300])
% hold off
% set(gcf,'position',[10,10,1400,1000]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig3a1.png')
% 
% % UTP
% figure(2)
% plot(t2,y2,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% title('UTP')
% hold on
% % h = vline2(250,'r--');
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1300])
% hold off
% set(gcf,'position',[10,10,1400,1000]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig3c1.png')
% % 
% % % ADP
figure(3)
plot(t3,y3,'k','LineWidth',4)
xlabel('time (s)')
ylabel({'Fura-2 ratio'})
title('ADP')
hold on
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',6)
ax.FontSize=70;
box off
ylim([0.9 1.46])
yticks([0.9 1 1.1 1.2 1.3 1.4])
xlim([0 1750])
hold off
set(gcf,'position',[10,10,1400,1000]) %[xpos, ypos, Width, Height]
saveas(gcf,'Fig3b1.png')
% 
% % ADP+BIM
% figure(4)
% plot(t4,y4,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% % xlim([0 1300])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig5a','epsc')
% 
% % UTP+BIM
% figure(5)
% plot(t5,y5,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% xline(1070,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1700])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig5b','epsc')
% 
% % ADP (PKC-DR)
% figure(6)
% plot(t6,y6,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1700])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig4a','epsc')
% 
% % figure(7)
% plot(t7,y7,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1850])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig4b','epsc')
% 
% % UTP (PKC-DR)
% figure(8)
% plot(t8,y8,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1800])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig4c','epsc')
% 
% % figure(9)
% plot(t9,y9,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1800])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig4d','epsc')
% 
% % ADP (+-Ca)
% figure(10)
% plot(t10,y10,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1210])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig6a','epsc')
% 
% % figure(11)
% plot(t11,y11,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1450])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig6b','epsc')
% 
% % UTP (+- Ca)
% figure(12)
% plot(t12,y12,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1400])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig6c','epsc')
% 
% figure(13)
% plot(t13,y13,'k','LineWidth',4)
% xlabel('time (s)')
% ylabel({'Fura-2 ratio'})
% hold on
% % xline(910,'r--','LineWidth',6)
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% box off
% % xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1400])
% hold off
% set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig6d','epsc')
