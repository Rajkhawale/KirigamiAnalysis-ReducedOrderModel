function InputData = Str_and_Mesh_Gen(InputData)

    % Compute distances
    Rt = 2 / (InputData.Nt * (InputData.LcRt+1)); % Distances between the cuts in the transverse direciton
    Lc = Rt * InputData.LcRt; % Cut length
    
    syms x y
    eqn1 = y + 2*InputData.Nl*x == 2;
    eqn2 = InputData.D2Rl * y == x;
    [A,B] = equationsToMatrix([eqn1, eqn2], [x, y]);
    X = linsolve(A,B);
    Rl = double(X(1)); % Distances between the cuts in the longitudinal direciton
    D1 = double(X(2)); % Cut height
    
    cell_width_1 = 3*(InputData.Geo_len/(InputData.Nl*2+1));
    cell_width_other = 2*(InputData.Geo_len/(InputData.Nl*2+1));
    cell_length = InputData.Geo_len/InputData.Nt;
    cell_margin = (InputData.Geo_wd - InputData.Geo_len)/2;% cell starts above certain y distance.
    
    % Outer rectangle (cell boundary)
    boundary = [3; 4; 0; InputData.Geo_len; InputData.Geo_len; 0; 0; 0; InputData.Geo_wd; InputData.Geo_wd];
    
    % Initialize gd
    gd = boundary;
    ns = {'B'};
    sf_parts = {'B'};
    part_count = 0;
    cut_polygons = {};    
    
    %% For first cell (bottom left corner)
    % For central downward facing cut
    
    for j = 1:InputData.Nl
        for i = 1:InputData.Nt
            part_count = part_count + 1;
            [xp, yp] = generate_thick_curved_cut(Rt/2+cell_length*(i-1), (cell_width_1-D1-Rl)+cell_width_other*(j-1), Lc, D1, InputData, false, cell_margin, 'none'); % For central downward cut
            gd = append_cut(gd, xp, yp);
            cut_polygons{end+1} = [xp(:), yp(:)];
            ns{end+1} = sprintf('P%d', part_count);
            sf_parts{end+1} = ['-P' num2str(part_count)];
            
            if j==1
                % For 4 half upward facing cuts
                corner_positions = [cell_length*(i-1), cell_width_1;
                                    cell_length+cell_length*(i-1), cell_width_1;
                                    cell_length*(i-1), cell_width_1-2*Rl;
                                    cell_length+cell_length*(i-1), cell_width_1-2*Rl];
                for k = 1:4
                    part_count = part_count + 1;
                    [xp, yp] = generate_thick_curved_cut(corner_positions(k,1), corner_positions(k,2), Lc, D1, InputData, true, cell_margin, k);
                    gd = append_cut(gd, xp, yp);
                    cut_polygons{end+1} = [xp(:), yp(:)];
                    ns{end+1} = sprintf('P%d', part_count);
                    sf_parts{end+1} = ['-P' num2str(part_count)];
                end
            else
                % For 2 half upward facing cuts
                corner_positions = [cell_length*(i-1), cell_width_1+cell_width_other*(j-1);
                                    cell_length+cell_length*(i-1), cell_width_1+cell_width_other*(j-1)];
                for m = 1:2
                    part_count = part_count + 1;
                    [xp, yp] = generate_thick_curved_cut(corner_positions(m,1), corner_positions(m,2), Lc, D1, InputData, true, cell_margin, m);
                    gd = append_cut(gd, xp, yp);
                    cut_polygons{end+1} = [xp(:), yp(:)];
                    ns{end+1} = sprintf('P%d', part_count);
                    sf_parts{end+1} = ['-P' num2str(part_count)];
                end
            end
        
        end
    end
    %% Build decsg input
    sf = strjoin(sf_parts, ' ');
    ns_matrix = char(ns{:})';
    
    % Build geometry
    [g, ~] = decsg(gd, sf, ns_matrix);
    model = createpde;
    geometryFromEdges(model, g);
    
    % Generate mesh
    Mesh = generateMesh(model, 'Hmax', InputData.mesh_size, 'GeometricOrder', 'linear');
    
    figure;
    pdemesh(model);
    axis equal;
    % title('Single Kirigami Cell Mesh');
    
    % Store mesh
    InputData.coordinates = [Mesh.Nodes' zeros(length(Mesh.Nodes'),1)];
    InputData.triangles = Mesh.Elements';
    fprintf('Mesh Generation completed!\n');

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
          'EdgeColor', 'none', ...
          'LineWidth', 2.2, ...
          'FaceAlpha', 0.9);
    hold on;
    % if size(V,2)==2
        X = [V(bE(:,1),1) V(bE(:,2),1)]';             % 2D
        Y = [V(bE(:,1),2) V(bE(:,2),2)]';
        plot(X, Y, 'k-', 'LineWidth', 2);
        axis equal off;

    %% Get Cut Face IDs
    % For each cut polygon
    num_cuts = length(cut_polygons);
    cut_node_ids = cell(num_cuts, 1);
    Original_cord_mat = cell(num_cuts, 1);
    % x_nodes = Mesh.Nodes(1,:)';
    % y_nodes = Mesh.Nodes(2,:)';
    A_mat = Mesh.Nodes';
    for i = 1:num_cuts
        poly = cut_polygons{i};
        [~, cut_node_ids{i}] = ismember(round(poly,6), round(A_mat,6), 'rows');

        % x_poly = poly(:,1); y_poly = poly(:,2);
        % in = inpolygon(x_nodes, y_nodes, x_poly, y_poly);
        % cut_node_ids{i} = find(in);

        % Formulate the matrix that has same X cord nodes on cuts. This is to calculate total cut open area
        Original_cord_mat{i} = [cut_node_ids{i}(1:length(cut_node_ids{i})/2), cut_polygons{i}(1:length(cut_node_ids{i})/2,1), cut_polygons{i}(1:length(cut_node_ids{i})/2,2), ...
                  cut_polygons{i}(length(cut_node_ids{i}):-1:length(cut_node_ids{i})/2+1,1), cut_polygons{i}(length(cut_node_ids{i}):-1:length(cut_node_ids{i})/2+1,2), ...
                  cut_node_ids{i}(length(cut_node_ids{i}):-1:length(cut_node_ids{i})/2+1), abs(cut_polygons{i}(1:length(cut_node_ids{i})/2,2) - cut_polygons{i}(length(cut_node_ids{i}):-1:length(cut_node_ids{i})/2+1,2))];

    end
    InputData.cut_node_ids = cut_node_ids;
    InputData.Original_cord_mat = Original_cord_mat;

    


    
end

%% --- Helper: Generate thick curved cut polygon ---
function [x_poly, y_poly] = generate_thick_curved_cut(x0, y0, Lc, D1, InputData, flipY, cell_margin, half_side)
    % x0 and y0 are the starting point for the curve.

    x_base = linspace(0, pi, InputData.curve_discretize);
    x_scaled = (x_base - min(x_base)) / (max(x_base)-min(x_base)) * Lc;

    y_base = D1 * sin(x_base);
    if flipY
        y_base = -y_base;
        x_scaled =x_scaled - Lc/2;
    end
    
    % Add cell margin
    y_base = y_base+cell_margin;

    % Inner curve
    x_in = x_scaled + x0;
    y_in = y_base + y0 -InputData.cut_thickness/2; % This ensure we are adding cut thickness in both side. 

    % Compute normal vectors
    % dx = [diff(x_in), 0];
    % dy = [diff(y_in), 0];
    % len = sqrt(dx.^2 + dy.^2);
    % len(len == 0) = 1;  % Prevent divide-by-zero
    % nx = -dy ./ len;
    % ny = dx ./ len;
    % % Outer curve
    % x_out = x_in + InputData.cut_thickness * nx;
    % y_out = y_in + InputData.cut_thickness * ny;

    % % Outer curve
    x_out = x_in;
    y_out = y_in + InputData.cut_thickness;

    if flipY && mod(half_side, 2) == 0 % If we need left side half (from cell point of view)
        x_in(ceil(length(x_in)/2)+1:end) = [];
        y_in(ceil(length(y_in)/2)+1:end) = [];
        x_out(ceil(length(x_out)/2)+1:end) = [];
        y_out(ceil(length(y_out)/2)+1:end) = [];
        % x_in = -Lc/2; x_out = -0.2 %Lc/2; 
    elseif flipY && mod(half_side, 2) == 1 % If we need right side half (from cell point of view)
        x_in(1:floor(length(x_in)/2)) = [];
        y_in(1:floor(length(y_in)/2)) = [];
        x_out(1:floor(length(x_out)/2)) = [];
        y_out(1:floor(length(y_out)/2)) = [];
        % x_in = -Lc/2; x_out = -0.2 %Lc/2;
    end

    % Combine into polygon
    x_poly = [x_out, fliplr(x_in)];
    y_poly = [y_out, fliplr(y_in)];

    % Remove consecutive duplicates
    keep = [true, any(diff([x_poly; y_poly],1,2),1)];
    x_poly = x_poly(keep);
    y_poly = y_poly(keep);
end


%% --- Helper: Append cut polygon to gd ---
function gd_new = append_cut(gd, x, y)
    nvert = numel(x);
    cut_shape = [2; nvert; x(:); y(:)];
    % Pad existing gd and new cut to same row count
    max_rows = max(size(gd,1), length(cut_shape));
    gd = [gd; zeros(max_rows - size(gd,1), size(gd,2))];
    cut_shape = [cut_shape; zeros(max_rows - length(cut_shape), 1)];
    gd_new = [gd, cut_shape];
end