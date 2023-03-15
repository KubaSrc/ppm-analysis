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
l_az = 30;
l_el = 20;
face_alpha = .4;
edge_alpha = .3;

%%

% Set up a figure for plotting
figure(1); clf;
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