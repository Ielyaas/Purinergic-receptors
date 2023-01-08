clear
close all
clc

%% LOAD FILES

load('ki1p2y1bfdiag.dat')
load('ki1p2y1bfset.dat')
load('ki1p2y2bfdiag.dat')
load('ki1p2y2bfset.dat')
load('ki1p2y2bfset2.dat')
load('p2y2bfdiagfinal.dat')


% P2Y1
PKC_1 = [0 0.05 0.1 0.15 0.2 0.25];
%          0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7...
%          0.75 0.8 0.85 0.9 0.95 1];
PLC_1 = [0.026 0.046 0.055 0.06 0.063 0.061];
%          0.064 0.062 0.057 0.052 0.047 0.043...
%          0.04 0.038 0.036 0.034 0.033 0.031 0.03 0.03 0.029 0.029];
PKC_1f = [0 0.05 0.1 0.11 0.115 0.07 0.06 0.04...
          0.03 0.02 0.01 0];
PLC_1f = [0.033 0.058 0.081 0.086 0.089 0.116 0.125...
          0.137 0.145 0.152 0.151 0.088];

% P2Y2
PKC_2 = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55];
PLC_2 = [0.022 0.015 0.01 0.0095 0.0087 0.008 0.0072 0.0068 0.0066 0.0064...
         0.0063 0.0061];
PKC_2f = [0.0395 0.04 0.05 0.1 0.15 0.2 0.25 0.26 0.27 0.26 0.25 0.2 0.15...
          0.1 0.05 0.04 0.0395];
PLC_2f = [0.024 0.023 0.018 0.011 0.01 0.01 0.011 0.011 0.011 0.013 0.013...
          0.015 0.016 0.02 0.024 0.025 0.024];
      


%% PLOTS

% figure(1)
% plot(ki1p2y1bfdiag(1:52,1),ki1p2y1bfdiag(1:52,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(53:381,1),ki1p2y1bfdiag(53:381,2),'k--','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(382:480,1),ki1p2y1bfdiag(382:480,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(481:598,1),ki1p2y1bfdiag(481:598,2),'b--','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(481:598,1),ki1p2y1bfdiag(481:598,3),'b--','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(599:699,1),ki1p2y1bfdiag(599:699,2),'g','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(599:699,1),ki1p2y1bfdiag(599:699,3),'g','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(700:end,1),ki1p2y1bfdiag(700:end,2),'b--','LineWidth',3)
% hold on
% plot(ki1p2y1bfdiag(700:end,1),ki1p2y1bfdiag(700:end,3),'b--','LineWidth',3)
% axis([0 0.12 0 0.7])
% xlabel('V_{PLC} (\muM/s)')
% ylabel('[Ca^{2+}]_i \muM')
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% % ax.YAxis.Visible = 'off';
% box off
% hold off
% % set(gca,'ytick',[],'Ycolor','w','box','off')
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'figure10','epsc')

% figure(2)
% plot(ki1p2y2bfdiag(1:79,1),ki1p2y2bfdiag(1:79,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(80:371,1),ki1p2y2bfdiag(80:371,2),'k--','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(372:434,1),ki1p2y2bfdiag(372:434,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(435:635,1),ki1p2y2bfdiag(435:635,2),'b--','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(435:635,1),ki1p2y2bfdiag(435:635,3),'b--','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(636:788,1),ki1p2y2bfdiag(636:788,2),'g','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(636:788,1),ki1p2y2bfdiag(636:788,3),'g','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(789:942,1),ki1p2y2bfdiag(789:942,2),'b--','LineWidth',3)
% hold on
% plot(ki1p2y2bfdiag(789:942,1),ki1p2y2bfdiag(789:942,3),'b--','LineWidth',3)
% axis([0 0.1 0 0.5])
% xlabel('V_{PLC} (\muM/s)')
% ylabel('[Ca^{2+}]_i \muM')
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% % ax.YAxis.Visible = 'off';
% box off
% hold off
% % set(gca,'ytick',[],'Ycolor','w','box','off')
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'figure10','epsc')

figure(3)
plot(ki1p2y1bfset(1:39636,1),ki1p2y1bfset(1:39636,2),'b','LineWidth',3)
hold on
% plot(ki1p2y1bfset(39642:40426,1),ki1p2y1bfset(39642:40426,2),'g','LineWidth',3)
hold on
plot(ki1p2y1bfset(40432:40617,1),ki1p2y1bfset(40432:40617,2),'k','LineWidth',3)
hold on
plot(ki1p2y1bfset(40618:41670,1),ki1p2y1bfset(40618:41670,2),'b','LineWidth',3)
hold on
% plot(ki1p2y1bfset(41676:41763,1),ki1p2y1bfset(41676:41763,2),'g','LineWidth',3)
hold on
plot(ki1p2y1bfset(41769:end,1),ki1p2y1bfset(41769:end,2),'k','LineWidth',3)
hold on
plot(PLC_1,PKC_1,'Color',rgb('ForestGreen'),'LineWidth',3)
hold on
plot(PLC_1f,PKC_1f,'Color',rgb('DodgerBlue'),'LineWidth',3)
hold on
yline(0.1,'k--','LineWidth',3)
hold on
plot(0.04,0.1,'o','MarkerSize',20,...
    'MarkerEdgeColor',rgb('OrangeRed'),...
    'MarkerFaceColor',rgb('OrangeRed'))
hold on
plot(0.04,0,'o','MarkerSize',20,...
    'MarkerEdgeColor',rgb('DarkMagenta'),...
    'MarkerFaceColor',rgb('DarkMagenta'))
% grid minor
axis([0 0.2 0 1])
xlabel('V_{PLC} (\muM/s)')
ylabel('V_{PKC} (\muM/s)')
ax=gca;
set(ax,'Linewidth',6)
ax.FontSize=70;
% ax.YAxis.Visible = 'off';
box off

axes('Position',[.4 .4 .49 .51])
plot(ki1p2y1bfdiag(1:52,1),ki1p2y1bfdiag(1:52,2),'k','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(53:381,1),ki1p2y1bfdiag(53:381,2),'k--','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(382:480,1),ki1p2y1bfdiag(382:480,2),'k','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(481:598,1),ki1p2y1bfdiag(481:598,2),'b--','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(481:598,1),ki1p2y1bfdiag(481:598,3),'b--','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(599:699,1),ki1p2y1bfdiag(599:699,2),'g','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(599:699,1),ki1p2y1bfdiag(599:699,3),'g','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(700:end,1),ki1p2y1bfdiag(700:end,2),'b--','LineWidth',3)
hold on
plot(ki1p2y1bfdiag(700:end,1),ki1p2y1bfdiag(700:end,3),'b--','LineWidth',3)
hold on
plot(0.006772,0.08718,'o','MarkerSize',15,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
hold on
plot(0.08303,0.14327,'o','MarkerSize',15,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
hold on
plot(0.09945,0.26926,'o','MarkerSize',15,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
hold on
text(0.0089,0.088,'HB','FontSize',30)
hold on
text(0.075,0.162,'HB','FontSize',30)
hold on
text(0.101,0.26926,'SNPO','FontSize',30)
axis([0 0.11 0 0.65])
xlabel('V_{PLC} (\muM/s)')
ylabel('[Ca^{2+}]_i \muM')
ax=gca;
set(ax,'Linewidth',5)
ax.FontSize=40;
% ax.YAxis.Visible = 'off';
box off
hold off
% set(gca,'ytick',[],'Ycolor','w','box','off')
set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
saveas(gcf,'P2Y1_bf','epsc')

% figure(4)
% plot(ki1p2y2bfset(1:587,1),ki1p2y2bfset(1:587,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(588:4951,1),ki1p2y2bfset(588:4951,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(4957:5227,1),ki1p2y2bfset(4957:5227,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(5233:5590,1),ki1p2y2bfset(5233:5590,2),'g','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(5591:9046,1),ki1p2y2bfset(5591:9046,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(9047:9635,1),ki1p2y2bfset(9047:9635,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(9641:10510,1),ki1p2y2bfset(9641:10510,2),'k','LineWidth',3)
% hold on
% plot(ki1p2y2bfset(10516:end,1),ki1p2y2bfset(10516:end,2),'g','LineWidth',3)
% grid minor
% axis([0 0.4 0 1])
% xlabel('V_{PLC} (\muM/s)')
% ylabel('V_{PKC} (\muM/s)')
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% % set(gca, 'XScale', 'log')
% % set(gca, 'YScale', 'log')
% % ax.YAxis.Visible = 'off';
% box off
% hold off
% % set(gca,'ytick',[],'Ycolor','w','box','off')
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% % saveas(gcf,'figure2','epsc')

% figure(5)
% plot(ki1p2y2bfset2(1:463,1),ki1p2y2bfset2(1:463,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset2(464:7492,1),ki1p2y2bfset2(464:7492,2),'b','LineWidth',3)
% hold on
% % plot(ki1p2y2bfset2(7498:7718,1),ki1p2y2bfset2(7498:7718,2),'k','LineWidth',3)
% hold on
% % plot(ki1p2y2bfset2(7724:8052,1),ki1p2y2bfset2(7724:8052,2),'g','LineWidth',3)
% hold on
% plot(ki1p2y2bfset2(8503:11069,1),ki1p2y2bfset2(8503:11069,2),'b','LineWidth',3)
% hold on
% plot(ki1p2y2bfset2(11070:11481,1),ki1p2y2bfset2(11070:11481,2),'b','LineWidth',3)
% hold on
% % plot(ki1p2y2bfset2(11487:12104,1),ki1p2y2bfset2(11487:12104,2),'k','LineWidth',3)
% hold on
% % plot(ki1p2y2bfset2(12110:end,1),ki1p2y2bfset2(12110:end,2),'g','LineWidth',3)
% hold on
% yline(0.1,'k--','LineWidth',3)
% hold on
% plot(0.015,0.1,'o','MarkerSize',20,...
%     'MarkerEdgeColor',rgb('OrangeRed'),...
%     'MarkerFaceColor',rgb('OrangeRed'))
% hold on
% plot(0.015,0,'o','MarkerSize',20,...
%     'MarkerEdgeColor',rgb('DarkMagenta'),...
%     'MarkerFaceColor',rgb('DarkMagenta'))
% plot(PLC_2,PKC_2,'Color',rgb('ForestGreen'),'LineWidth',3)
% plot(PLC_2f,PKC_2f,'Color',rgb('DodgerBlue'),'LineWidth',3)
% % grid minor
% axis([0 0.12 0 0.8])
% xlabel('V_{PLC} (\muM/s)')
% ylabel('V_{PKC} (\muM/s)')
% ax=gca;
% set(ax,'Linewidth',6)
% ax.FontSize=70;
% % set(gca, 'XScale', 'log')
% % set(gca, 'YScale', 'log')
% % ax.YAxis.Visible = 'off';
% box off
% hold off
% % set(gca,'ytick',[],'Ycolor','w','box','off')
% set(gcf,'position',[10,10,1750,1400]) %[xpos, ypos, Width, Height]
% 
% axes('Position',[.4 .42 .48 .49])
% plot(p2y2bfdiagfinal(1:77,1),p2y2bfdiagfinal(1:77,2),'k','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(78:369,1),p2y2bfdiagfinal(78:369,2),'k--','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(370:440,1),p2y2bfdiagfinal(370:440,2),'k','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(441:637,1),p2y2bfdiagfinal(441:637,2),'b--','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(441:637,1),p2y2bfdiagfinal(441:637,3),'b--','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(638:783,1),p2y2bfdiagfinal(638:783,2),'g','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(638:783,1),p2y2bfdiagfinal(638:783,3),'g','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(784:end,1),p2y2bfdiagfinal(784:end,2),'b--','LineWidth',3)
% hold on
% plot(p2y2bfdiagfinal(784:end,1),p2y2bfdiagfinal(784:end,3),'b--','LineWidth',3)
% hold on
% plot(0.007453,0.089096,'o','MarkerSize',15,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% plot(0.04946,0.14324,'o','MarkerSize',15,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0]+0.9)
% hold on
% text(0.0085,0.09,'HB','FontSize',30)
% hold on
% text(0.045,0.157,'HB','FontSize',30)
% axis([0 0.07 0 0.45])
% xlabel('V_{PLC} (\muM/s)')
% ylabel('[Ca^{2+}]_i \muM')
% ax=gca;
% set(ax,'Linewidth',5)
% ax.FontSize=40;
% % ax.YAxis.Visible = 'off';
% box off
% hold off
% % saveas(gcf,'P2Y2_bf','epsc')

