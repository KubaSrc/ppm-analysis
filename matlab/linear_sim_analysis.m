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
plot(100*T.Var1,100*T.Var2,'-','LineWidth',3,'DisplayName','Simulation')
disp("Sk: ")
disp(T.Var2(10)./T.Var1(10))


%  Lc = 360
E_out = load('./results/kinematic/med-ppm/LC-360_Eout.mat').D;
E_in_RMSE = readtable('./results/kinematic/med-ppm/LC-254','Sheet','Lengths').E;
E_out_RMSE = sqrt(sum(E_out.^2)./length(E_out));
Lc = readtable('./results/kinematic/med-ppm/LC-360.xlsx','Sheet','Lengths').Lc;
s = scatter(100*E_in_RMSE./Lc,100*E_out_RMSE./Lc,'r','filled','lineWidth',3,'DisplayName','$L_c$ = 360 mm');
s.SizeData = 150;
s.Marker = "diamond";

%  Lc = 254
E_out = load('./results/kinematic/med-ppm/LC-254_Eout.mat').D;
E_in_RMSE = readtable('./results/kinematic/med-ppm/LC-254','Sheet','Lengths').E;
E_out_RMSE = sqrt(sum(E_out.^2)./length(E_out));
Lc = readtable('./results/kinematic/med-ppm/LC-254.xlsx','Sheet','Lengths').Lc;
s = scatter(100*E_in_RMSE./Lc,100*E_out_RMSE./Lc,'r','filled','lineWidth',3,'DisplayName','$L_c$ = 254 mm');
s.SizeData = 150;
s.Marker = "^";

%  Lc = 206
E_out = load('./results/kinematic/med-ppm/LC-206_Eout.mat').D;
E_in_RMSE = readtable('./results/kinematic/med-ppm/LC-206','Sheet','Lengths').E;
E_out_RMSE = sqrt(sum(E_out.^2)./length(E_out));
Lc = readtable('./results/kinematic/med-ppm/LC-206.xlsx','Sheet','Lengths').Lc;
s = scatter(100*E_in_RMSE./Lc,100*E_out_RMSE./Lc,'r','filled','lineWidth',3,'DisplayName','$L_c$ = 206 mm');
s.SizeData = 150;
s.Marker = "o";


% Format 
xlabel("$\epsilon_{in}$ \%","FontSize",38,'interpreter','latex')
ylabel("$\epsilon_{out}$ \%","FontSize",38,'interpreter','latex')
set(gcf,'color','w');
set(gcf,'position',[0,0,800,600])
lg  = legend('interpreter','latex');
lg.FontSize = 22;
set(lg,'Box','off')


xlim([0,0.1])