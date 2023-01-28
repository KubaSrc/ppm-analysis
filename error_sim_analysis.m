close all; clear all; clc

T = readtable('./data/error_sim_2023_01_27-07_04_43');

% Reduce volume to surface of optimal (H/D)*. Figure out what the optimal Vs = (H/Lc)

XY = [T.G,T.HD];
[XY_unique,ia,ic] = unique(XY,'rows');

G_unique = T.G(ia);
HD_unique = T.HD(ia);

Sk_star = inf(size(HD_unique,1),1);
H_star = NaN(size(HD_unique,1),1);

for idx = 1:size(ic)
    if T.E_mean(idx) < Sk_star(ic(idx))
        G_unique(ic(idx)) = T.G(idx);
        HD_unique(ic(idx)) = T.HD(idx);
        Sk_star(ic(idx)) = T.E_mean(idx);
        H_star(ic(idx)) = T.H(idx);
    end
end

Sk_star(isinf(Sk_star)) = NaN;

%%
disp("H optimal");
H_star(Sk_star == min(Sk_star))

disp("G optimal");
G_unique(Sk_star == min(Sk_star))

disp("HD optimal");
HD_unique(Sk_star == min(Sk_star))

disp("min Sk*");
disp(min(Sk_star))

figure()
plot(T.H(T.HD==0.5),T.E_mean(T.H(T.HD==0.5)))

%% Plotting parameters
ms = 100;
label_size = 38;
color_bar_size = 32;
legend_size = 18;
title_size = 36;
axis_size = 26;
font_type = 'arial';
alpha = 0.05;
AZ = 50;
EL = 25;
fig  = figure(1); clf; hold on;
fig.Position  = [100 100 1000 700];

% Reshape data for surface plot
mesh_HD = reshape(HD_unique,90,90);
mesh_G = reshape(G_unique,90,90);
mesh_K_star = reshape(Sk_star,90,90);

% Plot surface, contour, and reach constrained space
surf(mesh_HD,mesh_G,log(mesh_K_star)-5,log(mesh_K_star),'FaceColor','interp','EdgeColor','none')
[C,h] = contour(mesh_HD,mesh_G,log(mesh_K_star),[-3,-2,-1,0,1,2,3],'k','LineWidth',1.5,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black','LabelSpacing',144)
y_nan = [10,50,100,150,200,250,300,350]; x_nan = [0.09,1,2,3,4]; [X_nan, Y_nan] = meshgrid(x_nan,y_nan); Z_nan = repmat([-100],[8,5]);
surf(X_nan,Y_nan,Z_nan,'lineWidth',.5,'FaceColor','#696969','FaceAlpha',0.15);

% Just more matlab formatting for plot
set(gca,'fontsize',axis_size);
xlabel("$\frac{H}{D}$",'fontSize',label_size,'interpreter','latex')
ylabel("$\gamma$",'fontSize',label_size,'interpreter','latex')
zlabel("$log(\bar{S_k})$",'fontSize',label_size,'interpreter','latex')
view([AZ EL])
colormap(flipud(brewermap([],'Spectral')));
cb = colorbar;
clim([-2.75,2])
set(cb,'FontSize',color_bar_size)
cb.Label.Interpreter = 'latex';
cb.Label.String = 'arg min $log(Sk)$';
view(0,90)
xlim([0,4])

% lg = legend('interpreter','latex');
% lg.FontSize = 24;
% set(lg,'Box','off')

% 
% view(0,75)
% shading interp
% lightangle(-45,30)
% h.FaceLighting = 'gouraud';
% h.AmbientStrength = 0.3;
% h.DiffuseStrength = 0.8;
% h.SpecularStrength = 0.9;
% h.SpecularExponent = 25;
% h.BackFaceLighting = 'unlit';