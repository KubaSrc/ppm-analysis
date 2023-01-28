close all; clear all; clc

fontname = 'Helvetica';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

T = readtable('./ks_linear_sim.csv');
T = T(~isnan(T.Var2),:);


figure(5); clf; hold on;
ax = gca;
ax.FontSize = 24;

% Plot simulation data
plot(T.Var1,T.Var2,'-','LineWidth',3,'DisplayName','Simulation')
disp("Sk: ")
disp(T.Var2(10)./T.Var1(10))

E_out_med = load('./results/kinematic/med-ppm/error_out.mat').D;
E_in_med = readtable('./results/kinematic/med-ppm/lengths.xlsx').R;
E_out_RMSE = sqrt(sum(E_out_med.^2)./length(E_out_med));
E_in_RMSE = sqrt(sum(E_in_med.^2)./length(E_in_med));
Lc = 245.4617;

scatter(E_in_RMSE./Lc,E_out_RMSE./Lc,100,'r^','filled','lineWidth',3,'DisplayName','$L_c$ = 245mm'),

% Format 
xlabel("$\frac{\epsilon_{in}}{L_c}$","FontSize",56,'interpreter','latex')
ylabel("$\frac{\epsilon_{out}}{L_c}$","FontSize",56,'interpreter','latex')
set(gcf,'color','w');
set(gcf,'position',[0,0,800,600])
lg = legend('interpreter','latex');
lg.FontSize = 22;
set(lg,'Box','off')