clear all; close all; clc

map = brewermap(3,'Set1');

figure(1); hold on;

scatter(7.5e-5,0,200,map(1,:),'filled');
scatter(0.31,0,200,map(2,:), 'filled');
scatter(2.1,0,200,map(3,:),'filled');
xlim([9.99e-6,10])
ylim([0,1])
ax = gca;
ax.YAxis.Visible = 'off';
set(ax, 'XScale', 'log')
ax.FontSize = 20; 
set(gca, 'FontName', 'Times')
set(gca,'LineWidth',2);
set(gcf,'position',[0,0,1400,1000])