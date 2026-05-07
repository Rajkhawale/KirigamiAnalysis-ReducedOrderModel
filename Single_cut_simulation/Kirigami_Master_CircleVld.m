% Reset workspace
clearvars; close all;
tic
% Add required paths to analysis functions and fold pattern functions
addpath('Functions')

%% Kirigami Master File - For Generalized geometry and cut.

% --- INPUT PARAMETERS ---
InputData.R = 3; % m % Radius of circular sheet (inches)
InputData.cut_width = 2.58;% m
InputData.cut_height = 0.1;% m
InputData.thickness_sheet = 0.003; % inch - For mylar
InputData.density = 0.05; %1420; % 
InputData.elasticModulus = 760000;%  2.5e+09;%  %N/m2  - For mylar 
InputData.numberIncrements = 300;
InputData.maxIterations = 50;
InputData.loadMagnitude = 0.08; % m (stretching)
InputData.restAngleRad = 0;
InputData.loadType = 'Displacement';
InputData.foldStep = false;
InputData.plotNodes = 'no';
InputData.plotEnergy = 'no';
InputData.plotDeformedShape = 'yes';
InputData.plotReferenceShape = 'yes';
InputData.plotIncrement ='end';

%% Geometry and Mesh generation
% Circle center
cx = 0; cy = 0;

% Outer boundary: Circle in decsg format [1; center_x; center_y; radius]
outer_circle = [1; cx; cy; InputData.R];

% Inner rectangle cut
InputData.w = InputData.cut_width; InputData.h = InputData.cut_height;
x1 = cx - InputData.w/2; x2 = cx + InputData.w/2; y1 = cy - InputData.h/2; y2 = cy + InputData.h/2;
inner_rect = [3; 4; x1; x2; x2; x1; y1; y1; y2; y2];

% Make both columns same length WITHOUT inserting NaNs into `gd`
pad_to = max(length(outer_circle), length(inner_rect));
outer_circle = [outer_circle; zeros(pad_to - length(outer_circle),1)];
inner_rect = [inner_rect; zeros(pad_to - length(inner_rect),1)];

% Construct the geometry matrix
gd = [outer_circle, inner_rect];
sf = 'C1 - R1';  % Circle minus rectangle
ns = char('C1','R1'); ns = ns';

% Decompose the geometry
[g, bt] = decsg(gd, sf, ns);

% --- MESH GENERATION ---
model = createpde;
geometryFromEdges(model, g);

% Check edge labels
% figure; pdegplot(model, 'EdgeLabels', 'on'); axis equal; title('Edge Labels');

% Mesh with refinement (adjust edge number as needed)
Mesh = generateMesh(model, 'Hmax', 0.25, 'Hedge', {3,0.08}, 'Hedge', {4,0.08}, 'GeometricOrder', 'linear');  % You can add 'Hedge' after checking edge number
% Mesh = generateMesh(model, 'Hmax', 0.23, 'GeometricOrder', 'linear'); 

% % Plot mesh
figure;
pdemesh(model);
axis equal;
title('Triangular Mesh: Circular Sheet with Central Rectangular Cut');
 
% Store mesh
InputData.coordinates = [Mesh.Nodes' zeros(length(Mesh.Nodes'),1)];
InputData.triangles = Mesh.Elements';

% Figure: Draw undeformed structre
V   = InputData.coordinates;        % [x y] or [x y z]
tri = InputData.triangles;          % Mx3

% ----- find boundary edges -----
E = [tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])];   % all edges
E = sort(E,2);                                     % undirected
[uniqE,~,ic] = unique(E,'rows');
cnt = accumarray(ic,1);                            % how many times each edge appears
bE = uniqE(cnt==1,:);       

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);  % Large high-res figure
% Define face color from hex
hex = '#F7E444';  % Light yellow color F0F0F0 -> Light grey  F7E444 --> Yellow
rgbColor = sscanf(hex(2:end), '%2x%2x%2x', [1 3]) / 255;
faceColors = repmat(rgbColor, size(InputData.triangles, 1), 1);

patch('Faces', InputData.triangles, 'Vertices', InputData.coordinates, ...
      'FaceVertexCData', faceColors, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'k', ...
      'LineWidth', 1.2, ...
      'FaceAlpha', 0.9);
hold on;
% if size(V,2)==2
    X = [V(bE(:,1),1) V(bE(:,2),1)]';             % 2D
    Y = [V(bE(:,1),2) V(bE(:,2),2)]';
    if size(V,2) == 3
        Z = [V(bE(:,1),3) V(bE(:,2),3)]';
    else
        Z = zeros(size(X)); % fallback if only 2D
    end
plot3(X, Y, Z, 'k-', 'LineWidth', 2);
axis equal off;


%% Run Analysis
% Preprocess - Structure gen (elements, bars), bonudary conditions, and stiffness formulation
[PreprocessData, InputData] = Pre_Process_kirigami(InputData);

% Perform analysis (Solver)
[PostprocessData] = Path_Analysis_Kirigami(InputData, PreprocessData);
% 
% % Postprocess output data
PostprocessData = Post_Process_Kirigami(InputData, PreprocessData, PostprocessData); 
% 
toc
% % Plot results
Plot_Results_Kirigami(InputData, PreprocessData, PostprocessData);

%% Plot Force-Displacement

% Abaqus_RF_data = load('Abaqus_RF_data.mat');
% Abaqus_Disp_data = load('Abaqus_Disp_data.mat');

only_y = PostprocessData.Fhis;
reac_force = sum(abs(only_y), 1); % Assuming Fhis is a matrix
% disp = linspace(0, 2*InputData.loadMagnitude, length(reac_force)); % Adjust size
disp = PostprocessData.Uhis(3*InputData.nodes_top(1)-1,:) *2;
strain = disp/(InputData.R*2);
f_c = (reac_force*InputData.R)/(InputData.elasticModulus * (InputData.thickness_sheet^3));

% Plot reaction force vs displacement
figure;
loglog(strain, f_c(f_c~= 0), '-', 'LineWidth', 4, 'MarkerSize', 8); % Use '-o' for markers to verify data points
% plot(disp/6, reac_force(reac_force~= 0), '-', 'LineWidth', 4, 'MarkerSize', 8); % Use '-o' for markers to verify data points

hold on
% plot(Abaqus_Disp_data.Abaqus_Disp_data, abs(Abaqus_RF_data.Abaqus_RF_data), '--', 'LineWidth', 4, 'MarkerSize', 8); 
xlabel('Strain', FontSize=19);
ylabel('Reaction Force (lbf)', FontSize=19);
% title('Reaction Force vs Displacement', FontSize=20);
ax = gca; % Get current axes
ax.FontSize = 16;
grid on;
% legend('Bar & Hinge Model', 'FE Model', FontSize=19);
% legend('show', 'Location', 'northwest'); % Other options: northeast, best, etc.

save("RF_plot_25__","strain", "reac_force")

%% PLOT BOUNDARY CONDITIONS 

% figure;
% n_all = scatter(InputData.coordinates(:,1),InputData.coordinates(:,2), 20, 'k');
% hold on
% n_Yfix = scatter(InputData.coordinates([InputData.nodes_bottom;InputData.nodes_top],1),InputData.coordinates([InputData.nodes_bottom;InputData.nodes_top],2), 90, 'b', '*', 'DisplayName', 'Y Fixed');
% hold on
% n_Xfix = scatter(InputData.coordinates([InputData.nodes_left;InputData.nodes_right],1),InputData.coordinates([InputData.nodes_left;InputData.nodes_right],2), 90, 'r', '*', 'DisplayName', 'X Fixed');
% hold on
% n_Zfix = scatter(InputData.coordinates([InputData.Zfix_node1; InputData.Zfix_node2],1),InputData.coordinates([InputData.Zfix_node1; InputData.Zfix_node2],2), 90, 'm','*', 'DisplayName', 'Z Fixed');
% legend([n_all, n_Yfix, n_Xfix, n_Zfix], {'All Nodes', 'Y Fixed', 'X Fixed', 'Z Fixed'}, FontSize=15); 

%% PLOT - Good resolution

figure('Color', 'w', 'Position', [100 100 1200 900]);  % White background and high resolution
filename = 'BuildingFacade.gif';
tri = InputData.triangles; % Triangular connectivity

% Fixed axis limits for consistent animation
all_coords = InputData.coordinates(:,1:3);
max_disp = max(vecnorm(PostprocessData.Uhis));  % Estimate max displacement
xlim_vals = [min(all_coords(:,1)), max(all_coords(:,1))];
ylim_vals = [min(all_coords(:,2)), max(all_coords(:,2))];
zlim_vals = [min(all_coords(:,3)), max(all_coords(:,3)) + max_disp];

for i = 1:size(PostprocessData.Uhis, 2)
    Ui = reshape(PostprocessData.Uhis(:,i), 3, [])';  % Deformation at increment
    deformed_coords = InputData.coordinates(:,1:3) + Ui;

    trisurf(tri, ...
        deformed_coords(:,1), ...
        deformed_coords(:,2), ...
        deformed_coords(:,3), ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'k', 'FaceAlpha', 0.9, 'LineWidth', 0.5);

    set(gca, 'LineWidth', 0.5);
    set(gcf, 'Renderer', 'opengl');   % Enables anti-aliasing
    title(sprintf('Increment: %d', i), 'FontSize', 14);
    axis equal;
    %xlim(xlim_vals); ylim(ylim_vals); %zlim(zlim_vals);
    view(3);
    shading interp;
    material dull;
    lighting gouraud;
    camlight('headlight');

    drawnow;

    % Save high-res frame and read it
    exportgraphics(gca, 'frame.png', 'Resolution', 100);
    img = imread('frame.png');
    [imind, cm] = rgb2ind(img, 256);

    % Write to GIF
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

% Optional: delete the temporary frame image
delete('frame.png');
