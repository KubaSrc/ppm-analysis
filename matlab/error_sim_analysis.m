close all; clear all; clc

T = readtable('./data/error_sim_2023_02_02-08_11_13');

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

%% Plot a histogram over full design space

fig  = figure(2); clf; hold on;
fig.Position  = [100 100 1000 700];

map = brewermap(3,'Set1');
color_idx = 1;

histogram(T.E_mean,926,'facecolor',"#A9A9A9",'EdgeColor','none')
histogram(T.E_mean,926,'EdgeColor',"#707070",'linewidth',2,'DisplayStyle','stairs')
histogram(T.E_mean,[0.06,0.075],'facecolor',map(color_idx,:),'EdgeColor','none')

xlim([0,2])

% Just more matlab formatting for plot
set(gca,'fontsize',axis_size);
xlabel("$S_k$",'fontSize',label_size,'interpreter','latex')
ylabel("Count",'fontSize',label_size,'interpreter','latex')


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
fig  = figure(2); clf; hold on;
fig.Position  = [100 100 1000 700];

% Reshape data for surface plot
mesh_HD = reshape(HD_unique,100,100);
mesh_G = reshape(G_unique,100,100);
mesh_K_star = reshape(Sk_star,100,100);

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
clim([min(log(Sk_star)),max(log(Sk_star))])
set(cb,'FontSize',color_bar_size)
cb.Label.Interpreter = 'latex';
cb.Label.String = 'min $log(Sk)$';
view(0,90)
xlim([0,3])

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
