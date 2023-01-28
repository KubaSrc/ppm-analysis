close all; clc; clear all;
set(0,'defaultfigurecolor',[1 1 1])

% Get data for theta-phi config
N = readtable("./nodes.csv");
N = table2array(N);
N = N.';
N_pairs = [[0;1],[1;2],[2;3],[2;4],[2;5],[3;6],[4;6],[5;6],[0;3],[0;4],[0;5],[3;5],[4;5]]+1;

% Get positions for home config
N_home = readtable("./nodes-home.csv");
N_home = table2array(N_home);
N_home = N_home.';

% Get positions for optimal config
N_optimal = readtable("./nodes-optimal.csv");
N_optimal = table2array(N_optimal);
N_optimal = 86.13.*N_optimal.';

% Get positions for theta config
N_theta = readtable("./nodes-theta.csv");
N_theta = table2array(N_theta);
N_theta = N_theta.';

% Get positions for phi config
N_phi = readtable("./nodes-phi.csv");
N_phi = table2array(N_phi);
N_phi = N_phi.';


% Plotting Variables
az = -15;
el = 20;
l_az = 20;
l_el = 20;
face_alpha = .4;
edge_alpha = .3;

%% Create plot for figure 1-a (reimann projection)

% Set up a figure for plotting
figure(1); clf;
hold on;
axis equal

% Variables for plotting surface projection
n = 30;
n_max = 30;
theta_max = deg2rad(35);
Lchar = N(1,7);
r = N(1,2);

% Sweep of max angles
phi = linspace(0,2*pi,n_max);
theta = linspace(0,theta_max,n_max);
[PHI,THETA] = meshgrid(phi,theta);


% Find input surface to project
X_in = r.*sin(THETA).*cos(PHI);
Y_in = r.*sin(THETA).*sin(PHI);
Z_in = r.*cos(THETA)+r;
[X,Y,Z] = sphere(n);
% surf(Z_in,Y_in,X_in,'FaceAlpha',face_alpha,"EdgeAlpha",edge_alpha,"FaceColor","#0072BD","LineWidth",1.5)
surf(.995.*Z.*r+r,.995.*Y.*r,.995.*X.*r,"FaceAlpha",0.1,'EdgeAlpha',0,'FaceColor',"#707070","LineWidth",1.5);

% % Find output surface
Zc = Lchar.*cot((pi-THETA)./2).*exp(i.*PHI);
Z_out = real(Zc);
Y_out = imag(Zc);
X_out = repmat([Lchar],n_max,n_max);
surf(X_out,Y_out,Z_out,'FaceAlpha',face_alpha,"EdgeAlpha",0.4,"FaceColor","#0072BD","LineWidth",1.75);


% Draw circles on main axes of the sphere
theta = linspace(0,pi,100);
theta_back = linspace(pi,2*pi,100);
r = 1.005.*r;
% Circle 1
X = r.*cos(theta_back)+r; Y = r.*sin(theta_back); Z = repmat([0],100);
plot3(X,Y,Z,'r','LineWidth',4)
X = r.*cos(theta)+r; Y = r.*sin(theta); Z = repmat([0],100);
plot3(X,Y,Z,'--r','LineWidth',4)
% Circle 2
X = r.*cos(theta)+r; Y = repmat([0],100); Z = r.*sin(theta);
plot3(X,Y,Z,'b','LineWidth',4)
X = r.*cos(theta_back)+r; Y = repmat([0],100); Z = r.*sin(theta_back);
plot3(X,Y,Z,'--b','LineWidth',4)
% Circle 3
X = repmat([0],100)+r; Y = r.*sin(theta_back); Z = r.*cos(theta_back);
plot3(X,Y,Z,'g','LineWidth',4)
X = repmat([0],100)+r; Y = r.*sin(theta); Z = r.*cos(theta);
plot3(X,Y,Z,'--g','LineWidth',4)


% Plot of line being projected
theta_trace = deg2rad(27.5);
phi_trace = deg2rad(35);
z_proj = Lchar.*cot((pi-theta_trace)/2).*exp(i.*phi_trace);
% Point on the sphere
zx = r.*cos(theta_trace); zy = r.*sin(theta_trace).*sin(phi_trace); zz = r.*sin(theta_trace)*cos(phi_trace);
% Projected point
zxp = Lchar; zyp = imag(z_proj); zzp = real(z_proj);
% Draw line of projection
plot3([0,zxp],[0,zyp],[0,zzp],'k--',LineWidth=2)
% Draw point on sphere
[X,Y,Z] = sphere();
rr = 1;
surf(zx+rr.*X+r,zy+rr.*Y,zz+rr.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
surf(zxp+rr.*X,zyp+rr.*Y,zzp+rr.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
% Points to label
surf(2.*r+rr.*X,rr.*Y,rr.*Z,'FaceColor',"black",'EdgeColor','none');
surf(r+rr.*X,r+rr.*Y,rr.*Z,'FaceColor',"black",'EdgeColor','none');
surf(r+rr.*X,-r+rr.*Y,rr.*Z,'FaceColor',"black",'EdgeColor','none');
surf(r+rr.*X,rr.*Y,r+rr.*Z,'FaceColor',"black",'EdgeColor','none');
surf(r+rr.*X,rr.*Y,-r+rr.*Z,'FaceColor',"black",'EdgeColor','none');
surf(rr.*X,rr.*Y,rr.*Z,'FaceColor',"black",'EdgeColor','none');

axis equal
set(gca,'visible','off')
view(az,el)
set(gcf,'position',[0,0,1400,1000])
lightangle(az,el)
lighting gouraud

camproj('perspective')
%% Create plot for figure 1-b (mechanism only)

% Set up a figure for plotting
figure(2); clf;
hold on; grid on


% Draw all of the links
for pair = N_pairs
    R = .5; S = 20;
    r1 = N_home(:,pair(1));
    r2 = N_home(:,pair(2));
    [X_cyl,Y_cyl,Z_cyl] = cylinder2P(R,S,r1.',r2.');
    surf(X_cyl,Z_cyl,Y_cyl,'FaceColor',"#A0A0A0",'EdgeColor','none')
end

% Draw all of the nodes
for node = N_home
    [X,Y,Z] = sphere();
    r = 1;
    surf(node(1)+r.*X,node(3)+r.*Y,node(2)+r.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
    
end


% Change view
axis equal
set(gca,'visible','off')
view(az,el)
set(gcf,'position',[0,0,1400,1000])
lightangle(l_az,l_el)
lighting gouraud

camproj('perspective')
camroll(-5)

%% Create plot for figures 1-c (theta-phi mechanism)

% Set up a figure for plotting
f = figure(3); clf;
hold on; grid on


% Draw all of the links
for pair = N_pairs
    R = .5; S = 20;
    r1 = N(:,pair(1));
    r2 = N(:,pair(2));
    [X_cyl,Y_cyl,Z_cyl] = cylinder2P(R,S,r1.',r2.');
    surf(X_cyl,Z_cyl,Y_cyl,'FaceColor',"#A0A0A0",'EdgeColor','none')
end

% Draw all of the nodes
for node = N
    [X,Y,Z] = sphere();
    r = 1;
    surf(node(1)+r.*X,node(3)+r.*Y,node(2)+r.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
    
end

% Variables for plotting surface projection
n = 20;
n_max = 25;
theta_max = deg2rad(35);
Lchar = N(1,7);
r = N(1,2);

% Sweep of max angles
phi = linspace(0,2*pi,n_max);
theta = linspace(0,theta_max,n_max);
[PHI,THETA] = meshgrid(phi,theta);


% Find input surface to project
X_in = r.*sin(THETA).*cos(PHI);
Y_in = r.*sin(THETA).*sin(PHI);
Z_in = r.*cos(THETA)+r;
[X,Y,Z] = sphere(n);
surf(Z_in,Y_in,X_in,'FaceAlpha',face_alpha,"EdgeAlpha",edge_alpha,"FaceColor","#0072BD","LineWidth",1.5)
surf(.995.*Z.*r+r,.995.*Y.*r,.995.*X.*r,"FaceAlpha",0.1,'EdgeAlpha',0,'FaceColor',"#707070","LineWidth",1.5);

% % Find output surface
Zc = Lchar.*cot((pi-THETA)./2).*exp(i.*PHI);
Z_out = real(Zc);
Y_out = imag(Zc);
X_out = repmat([Lchar],n_max,n_max);
surf(X_out,Y_out,Z_out,'FaceAlpha',face_alpha,"EdgeAlpha",edge_alpha,"FaceColor","#0072BD","LineWidth",1.5);

% Plot of line being projected
theta_trace = deg2rad(27.5);
phi_trace = deg2rad(35);
z_proj = Lchar.*cot((pi-theta_trace)/2).*exp(i.*phi_trace);
% Point on the sphere
zx = r.*cos(theta_trace); zy = r.*sin(theta_trace).*sin(phi_trace); zz = r.*sin(theta_trace)*cos(phi_trace);
% Projected point
zxp = Lchar; zyp = imag(z_proj); zzp = real(z_proj);
% Draw line of projection
plot3([0,zxp],[0,zyp],[0,zzp],'k--',LineWidth=2)


axis equal
set(gca,'visible','off')
view(az,el)
set(gcf,'position',[0,0,1400,1000])
lightangle(az,el)
lighting gouraud

camproj('perspective')

%% Labeling for Design Parametrization

% Set up a figure for plotting
figure(4); clf;
hold on; grid on

N_pairs = [[0;1],[1;2],[2;3],[2;4],[2;5],[3;6],[4;6],[5;6],[0;3],[0;4],[0;5],[3;5],[4;5]]+1;
l_alpha = 0.2;
N_alpha = [l_alpha,l_alpha,1,1,1,1,1,1,l_alpha,l_alpha,l_alpha,1,1];

% Draw all of the links
for pair_idx = 1:length(N_pairs)
    pair = N_pairs(:,pair_idx);
    R = .5; S = 20;
    r1 = N_home(:,pair(1));
    r2 = N_home(:,pair(2));
    [X_cyl,Y_cyl,Z_cyl] = cylinder2P(R,S,r1.',r2.');
    surf(X_cyl,Z_cyl,Y_cyl,'FaceColor',"#A0A0A0",'EdgeColor','none','FaceAlpha',N_alpha(pair_idx))
end

% Draw all of the nodes
for node = N_home
    [X,Y,Z] = sphere();
    r = 1;
    surf(node(1)+r.*X,node(3)+r.*Y,node(2)+r.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
    
end

% Draw lines for D
D_center = [(N_home(1,7) - N_home(1,3))/2 + N_home(1,3),0,0];
plot3([D_center(1),N_home(1,6)],[D_center(3),N_home(3,6)],[D_center(2),N_home(2,6)],'k','lineWidth',3)
plot3([D_center(1),N_home(1,5)],[D_center(3),N_home(3,5)],[D_center(2),N_home(2,5)],'k--','lineWidth',3)
plot3([D_center(1),N_home(1,4)],[D_center(3),N_home(3,4)],[D_center(2),N_home(2,4)],'k--','lineWidth',3)

G_theta = linspace(deg2rad(120),deg2rad(240),100);
G_Z = 6*cos(G_theta);
G_Y = 6*sin(G_theta);
G_X = repmat(D_center(1),100);
plot3(G_X,G_Y,G_Z,'k','lineWidth',3)


plot3([D_center(1),N_home(1,7)],[D_center(3),N_home(3,7)],[D_center(2),N_home(2,7)],'k','lineWidth',3)
% Change view
axis equal
set(gca,'visible','off')
view(32.5,5)
set(gcf,'position',[0,0,1400,1000])
lightangle(40,l_el)
lighting gouraud
camproj('perspective')

%% Create plot for figure 1-b (mechanism only)

% Set up a figure for plotting
figure(5); clf;
hold on; grid on


% Draw all of the links
for pair = N_pairs
    R = .5; S = 20;
    r1 = N_optimal(:,pair(1));
    r2 = N_optimal(:,pair(2));
    [X_cyl,Y_cyl,Z_cyl] = cylinder2P(R,S,r1.',r2.');
    surf(X_cyl,Z_cyl,Y_cyl,'FaceColor',"#A0A0A0",'EdgeColor','none')
end

% Draw all of the nodes
for node = N_optimal
    [X,Y,Z] = sphere();
    r = 1;
    surf(node(1)+r.*X,node(3)+r.*Y,node(2)+r.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
    
end


% Change view
axis equal
set(gca,'visible','off')
view(45,el)
set(gcf,'position',[0,0,1400,1000])
lightangle(l_az,l_el)
lighting gouraud

camproj('perspective')
camroll(-5)

%% Labeling for Design Parametrization

% Set up a figure for plotting
figure(6); clf;
hold on; grid on

N_pairs = [[0;1],[1;2],[2;3],[2;4],[2;5],[3;6],[4;6],[5;6],[0;3],[0;4],[0;5],[3;5],[4;5]]+1;
l_alpha = 0.2;
N_alpha = [l_alpha,l_alpha,1,1,1,1,1,1,l_alpha,l_alpha,l_alpha,1,1];

% Draw all of the links
for pair_idx = 1:length(N_pairs)
    pair = N_pairs(:,pair_idx);
    R = .5; S = 20;
    r1 = N_optimal(:,pair(1));
    r2 = N_optimal(:,pair(2));
    [X_cyl,Y_cyl,Z_cyl] = cylinder2P(R,S,r1.',r2.');
    surf(X_cyl,Z_cyl,Y_cyl,'FaceColor',"#A0A0A0",'EdgeColor','none','FaceAlpha',N_alpha(pair_idx))
end

% Draw all of the nodes
for node = N_optimal
    [X,Y,Z] = sphere();
    r = 1;
    surf(node(1)+r.*X,node(3)+r.*Y,node(2)+r.*Z,'FaceColor',"#FF0000",'EdgeColor','none');
    
end

% Draw lines for D
D_center = [(N_optimal(1,7) - N_optimal(1,3))/2 + N_optimal(1,3),0,0];
plot3([D_center(1),N_optimal(1,6)],[D_center(3),N_optimal(3,6)],[D_center(2),N_optimal(2,6)],'k','lineWidth',3)
plot3([D_center(1),N_optimal(1,5)],[D_center(3),N_optimal(3,5)],[D_center(2),N_optimal(2,5)],'k--','lineWidth',3)
plot3([D_center(1),N_optimal(1,4)],[D_center(3),N_optimal(3,4)],[D_center(2),N_optimal(2,4)],'k--','lineWidth',3)

G_theta = linspace(deg2rad(120+20),deg2rad(240-20),100);
G_Z = 6*cos(G_theta);
G_Y = 6*sin(G_theta);
G_X = repmat(D_center(1),100);
plot3(G_X,G_Y,G_Z,'k','lineWidth',3)
plot3([D_center(1),N_optimal(1,7)],[D_center(3),N_optimal(3,7)],[D_center(2),N_optimal(2,7)],'k','lineWidth',3)

% Change view
axis equal
set(gca,'visible','off')
view(32.5,5)
set(gcf,'position',[0,0,1400,1000])
lightangle(40,l_el)
lighting gouraud
camproj('perspective')
