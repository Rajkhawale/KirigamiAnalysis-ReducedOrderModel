% Reset workspace
clearvars; close all;
tic
% Add required paths to analysis functions and fold pattern functions
addpath('Functions')

%% Kirigami Master File - For Generalized geometry and cut.

% --- INPUT PARAMETERS ---
InputData.Geo_len = 2; % (m) Length of entire structure
InputData.Geo_wd = 2.4;% (m) Width of entire structure
InputData.thickness_sheet = 0.003; % (m) plastic [Original from Abaqus 0.002]
InputData.density = 1420;
InputData.elasticModulus = 2.5e+09; % N/m2
InputData.numberIncrements = 100;
InputData.maxIterations = 50;
InputData.restAngleRad = 0;
InputData.loadType = 'Displacement';
InputData.foldStep = false;
InputData.plotNodes = 'no';
InputData.plotEnergy = 'no';
InputData.plotDeformedShape = 'yes';
InputData.plotReferenceShape = 'yes';
InputData.plotIncrement ='end';

%% Cut configuration (INPUT)

% File = 192
Tnum = 192;
InputData.Nt = 5; % Nt - Number of cells in transverse direction
InputData.Nl = 10;% Nl - Number of cells in longitudinal direction
InputData.LcRt = 2;% (ratio) Lc - Length of cut; Rt - Dist. between cuts in transverse direction
InputData.D2Rl = 1;
% InputData.perturb = 0.09; % Initial displacement in z-direction. 
InputData.mesh_size = 0.032; 
InputData.curve_discretize = 13;
InputData.cut_thickness = 0.005;  
InputData.loadMagnitude = 0.4; % 0.44 in abaqus - Only in top direciton

%% Structure and mesh generation

InputData = Str_and_Mesh_Gen(InputData);

% Visualize undeformed mesh strcuture
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
[PreprocessData, InputData] = Pre_Process_BuildingFacade(InputData);

% Perform analysis (Solver)
[PostprocessData, InputData] = Path_Analysis_BuildingFacade(InputData, PreprocessData);

% Postprocess output data
PostprocessData = Post_Process_Kirigami(InputData, PreprocessData, PostprocessData); 
% 
toc

%% PLOT GIF - For Good resolution
% 
figure('Color', 'w', 'Position', [100 100 1200 900]);  % White background and high resolution
filename = 'Celular_20_5.gif';
tri = InputData.triangles; % Triangular connectivity

% Fixed axis limits for consistent animation
all_coords = InputData.coordinates(:,1:3);
max_disp = max(vecnorm(PostprocessData.Uhis));  % Estimate max displacement
xlim_vals = [min(all_coords(:,1)), max(all_coords(:,1))];
ylim_vals = [min(all_coords(:,2)), max(all_coords(:,2))];
zlim_vals = [min(all_coords(:,3)), max(all_coords(:,3)) + max_disp];

% Rotation matrix around Y-axis for 40 degrees
angle_deg = 40;
angle_rad = deg2rad(angle_deg);

for i = 1:size(PostprocessData.Uhis, 2)
    Ui = reshape(PostprocessData.Uhis(:,i), 3, [])';  % Deformation at increment
    deformed_coords = InputData.coordinates(:,1:3) + Ui;

    % Apply rotation
    rotated_coords = deformed_coords;% * R_y;

    trisurf(tri, ...
        rotated_coords(:,1), ...
        rotated_coords(:,2), ...
        rotated_coords(:,3), ...
        rotated_coords(:,2), ...
        'FaceColor', 'interp', ...
        'EdgeColor', 'k', 'FaceAlpha', 0.9, 'LineWidth', 0.5);

    set(gca, 'LineWidth', 0.5);
    set(gcf, 'Renderer', 'opengl');   % Enables anti-aliasing
    % title(sprintf('Increment: %d', i), 'FontSize', 14);
    axis equal;
    %xlim(xlim_vals); ylim(ylim_vals); %zlim(zlim_vals);s
    % view(1);
    shading interp;
    material dull;
    lighting gouraud;
    camlight('headlight');
    grid off; 
    axis off; 
    colormap("summer"); 
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
% 
% Optional: delete the temporary frame image
delete('frame.png');

%% Open area calculation - From top view (Top view area is considered for airflow calculations.)
% 
% A_CutOpen_incr = zeros(size(PostprocessData.Uhis, 2),2); % As we have a total of 20% stretching.
% A_CutOpen_incr(:,2) = round(0:size(PostprocessData.Uhis, 2)/(size(PostprocessData.Uhis, 2)-1):size(PostprocessData.Uhis, 2));
% A_CutOpen_incr(1,2) = 1;
% 
% for pp = 1:length(A_CutOpen_incr)
% 
%     Ui = reshape(PostprocessData.Uhis(:,A_CutOpen_incr(pp,2)), 3, [])';  % [X Y Z] displacements, size: num_nodes × 3 % Extract displacement at current increment
%     deformed_coords = InputData.coordinates(:,1:3) + Ui; % New positions = original + displacement
%     A_CutOpen = zeros(length(InputData.cut_node_ids),1);
% 
%     for i = 1:length(InputData.cut_node_ids)
%         matrix = InputData.Original_cord_mat{i};
%         Def_cord_mat = [matrix(:,1), deformed_coords(matrix(:,1),1:2), deformed_coords(matrix(:,6),1:2), matrix(:,6)];
%         trapezoid_area = zeros(size(matrix,1)-1,1);
% 
%         for jj = 1:size(matrix,1)-1 % This is to loop over each trapezoid in cut.
% 
%             x = [matrix(jj,2), matrix(jj+1,2), matrix(jj+1,4), matrix(jj,4)]; % Original (undeformed) X coords
%             y = [Def_cord_mat(jj,3), Def_cord_mat(jj+1,3), Def_cord_mat(jj+1,5)+matrix(jj+1,7), Def_cord_mat(jj,5)+matrix(jj,7)]; % Substracting initial cut area
% 
%             trapezoid_area(jj) = 0.5 * abs( x(1)*y(2) + x(2)*y(3) + x(3)*y(4) + x(4)*y(1) - y(1)*x(2) - y(2)*x(3) - y(3)*x(4) - y(4)*x(1) );
%         end
%         A_CutOpen(i) = sum(trapezoid_area);
%     end
%     A_CutOpen_incr(pp,1) = sum(A_CutOpen);
% end
% 

    %% CUT OPEN AREA CALCULATION - Different Method

A_CutOpen_incr = zeros(size(PostprocessData.Uhis, 2),2); % As we have a total of 20% stretching.
A_CutOpen_incr(:,2) = round(0:size(PostprocessData.Uhis, 2)/(size(PostprocessData.Uhis, 2)-1):size(PostprocessData.Uhis, 2));
A_CutOpen_incr(1,2) = 1;

for pp = 1:length(A_CutOpen_incr)

    Ui = reshape(PostprocessData.Uhis(:,A_CutOpen_incr(pp,2)), 3, [])';  % [X Y Z] displacements, size: num_nodes × 3 % Extract displacement at current increment
    deformed_coords = InputData.coordinates(:,1:3) + Ui; % New positions = original + displacement
    A_CutOpen = zeros(length(InputData.cut_node_ids),1);

    % Ui = reshape(PostprocessData.Uhis(:,A_CutOpen_incr(end,2)), 3, [])';  % [X Y Z] displacements, size: num_nodes × 3 % Extract displacement at current increment
    % deformed_coords = InputData.coordinates(:,1:3) + Ui; % New positions = original + displacement

        % Step 1: Extract all edges from triangles
    edges = [InputData.triangles(:,[1 2]);
             InputData.triangles(:,[2 3]);
             InputData.triangles(:,[3 1])];
    
    % Sort edges to make them direction-independent
    edges_sorted = sort(edges, 2);
    
    % Step 2: Identify boundary edges (edges appearing only once)
    [unique_edges, ~, ic] = unique(edges_sorted, 'rows');
    counts = accumarray(ic, 1);
    
    boundary_edges = unique_edges(counts == 1, :);
    
    % Step 3: Build graph of boundary edges
    G = graph(boundary_edges(:,1), boundary_edges(:,2));
    
    % Step 4: Find connected components (each = one loop)
    bins = conncomp(G);
    numLoops = max(bins);
    
    loops = cell(numLoops,1);
    
    for i = 1:numLoops
        loop_nodes = find(bins == i);
        
        % Subgraph (nodes are reindexed internally as 1:length(loop_nodes))
        subG = subgraph(G, loop_nodes);
        
        % Start from local node index = 1 (NOT loop_nodes(1))
        ordered_local = dfsearch(subG, 1);
        
        % Map back to global node indices
        loops{i} = loop_nodes(ordered_local);
    end
    
    % Step 5: Compute area of each loop (project to XY)
    areas = zeros(numLoops,1);
    nodes = deformed_coords;
    
    for i = 1:numLoops
        pts = nodes(loops{i}, 1:2);
        
        % Ensure proper polygon (auto-fix self-intersections)
        pg = polyshape(pts(:,1), pts(:,2), 'Simplify', true);
        
        areas(i) = area(pg);
    end
    
    % Step 6: Identify outer boundary (largest area)
    [~, idx_outer] = max(abs(areas));
    
    outer_boundary = loops{idx_outer};
    
    % Step 7: Remaining loops are cuts (holes)
    cut_loops = loops;
    cut_loops(idx_outer) = [];
    
    % Step 8: Compute total open area using UNION of all holes
    pg_total = polyshape();
    
    for i = 1:length(cut_loops)
        pts = nodes(cut_loops{i}, 1:2);
        
        pg = polyshape(pts(:,1), pts(:,2), 'Simplify', true);
        
        pg_total = union(pg_total, pg);
    end
    
    A_open = area(pg_total);
    A_CutOpen_incr(pp,1) = A_open;

end

%% PLOT OPEN AREA CALCULATION FOR VERIFICATION
pg_total = polyshape();

figure; hold on; axis equal;

% Plot all cut polygons
for i = 1:length(cut_loops)
    pts = nodes(cut_loops{i}, 1:2);
    
    pg = polyshape(pts(:,1), pts(:,2), 'Simplify', true);
    
    % Plot individual cut
    plot(pg, 'FaceColor', 'red', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    
    % Union
    pg_total = union(pg_total, pg);
end


pts_outer = nodes(outer_boundary, 1:2);
pg_outer = polyshape(pts_outer(:,1), pts_outer(:,2));

plot(pg_outer, 'FaceColor', 'none', 'EdgeColor', 'blue', 'LineWidth', 2);

plot(pg_total, 'FaceColor', 'green', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

legend('Cuts (individual)', 'Outer boundary', 'Union (final)');
title(['Total Open Area = ', num2str(area(pg_total))]);

% Total area after adding side outer area.
Area_final = area(pg_total) + ( (InputData.Geo_len*(InputData.Geo_wd+InputData.loadMagnitude)) - area(pg_outer));
