%% Experimental data for figures

clear
close all
clc

%% DATA

[~,sheet_name] = xlsfinfo('Data.xlsx');
for k=1:numel(sheet_name)
  [~,~,data{k}] = xlsread('Data.xlsx',sheet_name{k});
end

%% PLOTS

figure(1)
plot(data{1}{:,1},data{1}{:,2},'k','LineWidth',4)
xlabel('time (s)')
ylabel({'Fura-2 ratio'})
hold on
% h = vline2(250,'r--');
ax=gca;
set(ax,'Linewidth',6)
ax.FontSize=70;
box off
% xticks([0 200 400 600 800 1000 1200 1400])
% xlim([0 1400])
hold off
set(gcf,'position',[10,10,2000,1400]) %[xpos, ypos, Width, Height]
% saveas(gcf,'Fig3a','epsc')