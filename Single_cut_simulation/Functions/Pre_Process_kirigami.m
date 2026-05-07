

function [PreprocessData, InputData] = Pre_Process_kirigami(InputData)
    %% Add geometric perturbation
    % target_point =  [0,0.18];% [0.0209, 0.7791]; %[0,0.18];%       % Center of perturbation
    % max_disp = 0.1;                  % Max Z displacement at center
    % num_nodes = 5;                   % Number of nearest nodes to perturb
    % 
    % % Compute distances from target point
    % dists = vecnorm(InputData.coordinates(:, 1:2) - target_point, 2, 2);
    % 
    % % Find the 7 closest nodes
    % 
    % [~, sorted_idx] = sort(dists);
    % perturb_nodes = sorted_idx(1:num_nodes);
    % 
    % % Normalize distances to [0, 1] range for weighting
    % local_dists = dists(perturb_nodes);
    % weights = 1 - local_dists / max(local_dists);  % linear decay: center = 1, farthest = 0
    % 
    % % Apply weighted Z displacement
    % InputData.coordinates(perturb_nodes, 3) = max_disp * weights;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%% Z perturbation: max at cut edge and decreasing till top portion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % max_disp = 0.1;             % Maximum Z displacement
    % y_cut = 0.05;               % Y location of the top cut line
    % x_min = -1.29; x_max = 1.29;
    % y_max = max(InputData.coordinates(:,2));  % Topmost Y in the structure
    % 
    % % Select all nodes above the cut (including those on the cut line)
    % nodes_above_cut = find(InputData.coordinates(:,2) >= y_cut & ...
    %                        InputData.coordinates(:,1) >= x_min & ...
    %                        InputData.coordinates(:,1) <= x_max);
    % 
    % % Get Y-coordinates of selected nodes
    % y_vals = InputData.coordinates(nodes_above_cut, 2);
    % 
    % % Normalize Y from [y_cut → y_max] to [1 → 0]
    % weights = 1 - (y_vals - y_cut) / (y_max - y_cut);
    % weights = max(0, weights);  % ensure no negative weights
    % 
    % % Apply Z perturbation decreasing upward
    % InputData.coordinates(nodes_above_cut, 3) = max_disp * weights;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%% Z perturbation: max at st. line from cut to top point and gradually decreasing sideways %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Parameters
    max_disp = 0.05;             % Maximum Z displacement
    x_center = 0;               % X center of cut
    y_cut = 0.05;               % Y location of the top of the cut
    y_max = max(InputData.coordinates(:,2));  % Top of circular sheet

    % Get all nodes above the cut
    coords = InputData.coordinates(:, 1:2);
    nodes_above_cut = find(coords(:,2) >= y_cut);

    % Extract X and Y of those nodes
    x_vals = coords(nodes_above_cut, 1);
    y_vals = coords(nodes_above_cut, 2);

    % Normalize distances from cut center (for Gaussian-like decay)
    x_spread = 1.0;  % Controls lateral decay width
    y_spread = y_max - y_cut;  % Full vertical span

    % Distance from centerline (in x)
    dx = (x_vals - x_center) / x_spread;
    dy = (y_vals - y_cut) / y_spread;

    % 2D Gaussian weights: centered at (x=0, y=y_cut)
    weights = exp(-dx.^2 - dy.^2);

    % Apply Z perturbation
    InputData.coordinates(nodes_above_cut, 3) = max_disp * weights;
    
    
    
    
    %% Obtaining Bars
    % obtaining the bar elements
    InputData.bar_elems =  [InputData.triangles(:,1:2), InputData.triangles(:,2:3), InputData.triangles(:,[1,3])];
    InputData.bar_elems = [InputData.bar_elems(:,1:2); InputData.bar_elems(:,3:4); InputData.bar_elems(:,5:6)];
    InputData.bar_elems = unique(sort(InputData.bar_elems, 2), 'rows', 'stable');
    
    % Search for boundaries
    Edge = sort([InputData.triangles(:,1) InputData.triangles(:,2); InputData.triangles(:,2) InputData.triangles(:,3); InputData.triangles(:,3) InputData.triangles(:,1)],2);
    [u,~,n] = unique(Edge ,'rows');
    counts = accumarray(n(:), 1); 
    InputData.Bdry_bars = u(counts==1,:);
    
    
    
    %% Boundary conditions %%%%%%%% UPDATE %%%%%
    
    % [~, idx_sorted] = sort(InputData.coordinates(:,2));
    % nodes_bottom = idx_sorted(1:5); % Bottom-most nodes
    % nodes_top = idx_sorted(end-5+1:end); % Top-most nodes
    % 
    % [~, idx_sorted_x] = sort(InputData.coordinates(:,1));
    % nodes_left = idx_sorted_x(1:5); % Left-most nodes
    % nodes_right = idx_sorted_x(end-5+1:end); % Right-most nodes
    % 
    % % InputData.supports = [nodes_bottom, ones(length(nodes_bottom),1), ones(length(nodes_bottom),1), ones(length(nodes_bottom),1)];
    % % InputData.loads = [nodes_top, zeros(length(nodes_top),1), InputData.loadMagnitude*ones(length(nodes_top),1), zeros(length(nodes_top),1)];
    % %nodes_bottom, ones(length(nodes_bottom),1), ones(length(nodes_bottom),1), ones(length(nodes_bottom),1);
    % 
    % InputData.supports = [nodes_left, ones(length(nodes_left),1), ones(length(nodes_left),1), ones(length(nodes_left),1);
    %                         nodes_right, ones(length(nodes_right),1), ones(length(nodes_right),1), ones(length(nodes_right),1)];
    % 
    % % InputData.loads = [nodes_top, zeros(length(nodes_top),1), InputData.loadMagnitude*ones(length(nodes_top),1), zeros(length(nodes_top),1)];
    % 
    % 
    % % For stretching in both difection
    % % [~, idx_sorted_x] = sort(InputData.coordinates(:,1));
    % % nodes_left = idx_sorted_x(1:3); % Left-most nodes
    % % target_point1 = [-1.36,0]; target_point2 = [1.36,0]; target_point3 = [1.237, 0.443];
    % % dists1 = vecnorm(InputData.coordinates(:, 1:2) - target_point1, 2, 2);[~, load_node1] = min(dists1);
    % % dists2 = vecnorm(InputData.coordinates(:, 1:2) - target_point2, 2, 2);[~, load_node2] = min(dists2);
    % % dists3 = vecnorm(InputData.coordinates(:, 1:2) - target_point3, 2, 2);[~, load_node3] = min(dists3);
    % % 
    % % InputData.supports = [nodes_bottom, ones(length(nodes_bottom),1), zeros(length(nodes_bottom),1), zeros(length(nodes_bottom),1); % Fixing X bottom
    % %                        nodes_top, ones(length(nodes_top),1), zeros(length(nodes_top),1), zeros(length(nodes_top),1); % Fixing X top
    % %                     nodes_left, zeros(length(nodes_left),1), ones(length(nodes_left),1), zeros(length(nodes_left),1); % Fixing Y left
    % %                     nodes_right, zeros(length(nodes_right),1), ones(length(nodes_right),1), zeros(length(nodes_right),1); % Fixing Y right
    % %                     [load_node1; load_node2; load_node3], zeros(3,1), zeros(3,1), ones(3,1)]; % Fixing Z
    % % 
    % InputData.loads = [nodes_top, zeros(length(nodes_top),1), InputData.loadMagnitude*ones(length(nodes_top),1), zeros(length(nodes_top),1);
    %                     nodes_bottom, zeros(length(nodes_bottom),1), -InputData.loadMagnitude*ones(length(nodes_bottom),1), zeros(length(nodes_bottom),1)];
    
    %% Boundary conditions
    [~, idx_sorted] = sort(InputData.coordinates(:,2));
    InputData.nodes_bottom = idx_sorted(1:5); % Bottom-most nodes
    InputData.nodes_top = idx_sorted(end-5+1:end); % Top-most nodes
    
    [~, idx_sorted_x] = sort(InputData.coordinates(:,1));
    InputData.nodes_left = idx_sorted_x(1:5); % Left-most nodes
    InputData.nodes_right = idx_sorted_x(end-5+1:end); % Right-most nodes
    
    target_point1 = [-1.161,0.5]; target_point2 = [1.161,0.5];
    dists1 = vecnorm(InputData.coordinates(:, 1:2) - target_point1, 2, 2);[~, InputData.Zfix_node1] = min(dists1);
    dists2 = vecnorm(InputData.coordinates(:, 1:2) - target_point2, 2, 2);[~, InputData.Zfix_node2] = min(dists2);
    
    InputData.supports = [InputData.nodes_left, zeros(length(InputData.nodes_left),1), ones(length(InputData.nodes_left),1), zeros(length(InputData.nodes_left),1);
                            InputData.nodes_right, zeros(length(InputData.nodes_right),1), ones(length(InputData.nodes_right),1), zeros(length(InputData.nodes_right),1);
                            InputData.nodes_top, ones(length(InputData.nodes_top),1), zeros(length(InputData.nodes_top),1), zeros(length(InputData.nodes_top),1);
                            InputData.nodes_bottom, ones(length(InputData.nodes_bottom),1), zeros(length(InputData.nodes_bottom),1), zeros(length(InputData.nodes_bottom),1);
                            [InputData.Zfix_node1; InputData.Zfix_node2], zeros(2,1), zeros(2,1), ones(2,1)]; % Fixing Z
    
    InputData.loads = [InputData.nodes_top, zeros(length(InputData.nodes_top),1), InputData.loadMagnitude*ones(length(InputData.nodes_top),1), zeros(length(InputData.nodes_top),1);
                        InputData.nodes_bottom, zeros(length(InputData.nodes_bottom),1), -InputData.loadMagnitude*ones(length(InputData.nodes_bottom),1), zeros(length(InputData.nodes_bottom),1)];
    
    
    Ne = size(InputData.bar_elems,1); % Number of bars 
    Nn = size(InputData.coordinates,1); % Number of nodes.
    
    % Adding constraints in DOF form. rs -> FixedDOF 
    if size(InputData.supports,1) == 0
        rs = []; 
    else
        rs = [reshape([InputData.supports(:,1)*3-2,InputData.supports(:,1)*3-1,InputData.supports(:,1)*3]',[],1),...
              reshape(InputData.supports(:,2:4)',[],1)];
        rs(rs(:,2)==0,:)=[]; rs = rs(:,1);   
    end 
    PreprocessData.fixedDOFs = unique(rs);
    
    % Adding InputData.loads in DOF form
    if ~isempty(InputData.loads)
        FD = zeros(3*Nn,1);
        indp = InputData.loads(:,1);
        FD(3*indp-2) = InputData.loads(:,2); % X-load
        FD(3*indp-1) = InputData.loads(:,3); % Y-load
        FD(3*indp) = InputData.loads(:,4);   % Z-load
        DOF_load = FD;  % Puts the load back into the input options?
    end
    PreprocessData.DOF_load = DOF_load;
    
    %% Direction calculations -----> ??    
    
    Direction = [InputData.coordinates(InputData.bar_elems(:,2),1) - InputData.coordinates(InputData.bar_elems(:,1),1), ...  % containing the direction vector components [dx, dy, dz] for each bar.
         InputData.coordinates(InputData.bar_elems(:,2),2) - InputData.coordinates(InputData.bar_elems(:,1),2),...
         InputData.coordinates(InputData.bar_elems(:,2),3) - InputData.coordinates(InputData.bar_elems(:,1),3)];
    
    L_bar = sqrt(Direction(:,1).^2+Direction(:,2).^2+Direction(:,3).^2); % Length of each bar
    Direction = [Direction(:,1)./L_bar Direction(:,2)./L_bar Direction(:,3)./L_bar]; % Unit direction vector
    B = sparse(repmat((1:Ne)',1,6), ... %sparse matrix to relate bar forces to nodal displacements.
               [3*InputData.bar_elems(:,1)-2, 3*InputData.bar_elems(:,1)-1, 3*InputData.bar_elems(:,1), ...
                3*InputData.bar_elems(:,2)-2, 3*InputData.bar_elems(:,2)-1, 3*InputData.bar_elems(:,2)], ...
               [Direction, -Direction], Ne, 3*Nn);
    B = -B;
    PreprocessData.B = B;
    PreprocessData.L_bar = L_bar;
    
    %% Calculate bar areas
    % Shear coefficients
    a0 = 1.25; a1 = 0; a2 = 0.5;  %a0 = 0.45; a1 = 0; a2 = 0.4;
           
    % Effective Bar Areas
    barAreas = zeros(size(InputData.bar_elems,1),1);
    for i = 1:length(InputData.triangles) % Loop through panels
        % Find nodes on the ith panel
        for j = 1:3
            Index(j) = InputData.triangles(i,j);
            Point{j} = InputData.coordinates(Index(j),:);
        end; clear j
        % Calculate the area of the ith panel
        panelArea = 1/2*(norm(cross((Point{2}-Point{1}),(Point{3}-Point{1}))));
        % Calculate the length of the bars surrounding the ith panel
        for j = 1:3
            L(j) = norm(Point{j} - Point{j+1 - floor(j/3)*3});
        end; clear j
        Lengths = sort([L(1); L(2); L(3)]);
        AspectRatio(i) = Lengths(2)/Lengths(1); % As we are considering only alpha>1.
        % Calculate bar cross-sectional areas in the ith panel
        for j = 1:3
            if L(j) == Lengths(3)
                A(j) = (a2*AspectRatio(i)^2 + a1*AspectRatio(i) + a0)*panelArea*InputData.thickness_sheet/(L(j)); 
            else
                A(j) = panelArea*InputData.thickness_sheet/(L(j)); % factor value = 0.6
            end
        end; clear j
        % Assign panel bars area to final (global) bar_elems vector 
        % for j = 1:3
        %     endCondition = 0;
        %     k = 0;
        %     while endCondition == 0
        %         k = k + 1;
        %         find1 = ~isempty(find(InputData.bar_elems(k,:)==Index(j),1));
        %         find2 = ~isempty(find(InputData.bar_elems(k,:)==Index(j+1 - floor(j/3)*3),1));
        %         finder = find1 + find2;
        %         if finder == 2
        %             Bar_No(j) = k;
        %             endCondition = 1;
        %         end
        %     end
        % end; clear j

         for j = 1:3
        node1 = Index(j);
        node2 = Index(mod(j,3)+1);  % same as j+1 - floor(j/3)*3, cycles through 1→2→3→1
        
        % Check both possible node orders in bar_elems
        bars = InputData.bar_elems;
        matches = (bars(:,1)==node1 & bars(:,2)==node2) | ...
                  (bars(:,1)==node2 & bars(:,2)==node1);
        
        Bar_No(j) = find(matches, 1);  % Only first match, same as before
        end

        Panel_Bar_Areas = zeros(size(InputData.bar_elems,1),1);
        for j = 1:3
            Panel_Bar_Areas(Bar_No(j),1) = A(j);
        end
        barAreas = barAreas + Panel_Bar_Areas;
    end
    PreprocessData.barAreas = barAreas;
    PreprocessData.aspectRatio = AspectRatio;
    clear endCondition Bar_No find1 find2 finder Lengths k j Panel_Bar_Areas a0 a1 a2 A L
    
    %% Formulate bend bar matrix and Obtain bending stiffness array (Also initial angles)
    
    Comm = sparse(length(InputData.coordinates),size(InputData.triangles,1)); % This is to find triangular meshes that share two common nodes.
    for i=1:size(InputData.triangles,1) 
        Comm(InputData.triangles(i,:),i) = true;
    end
    Ge = Comm'*Comm;
    [mf, me] = find(triu(Ge==2)); % triangular meshes that share two common nodes
    Fold = [];
    Bend = [];
    for i=1:length(mf)
        [link,ia,ib] = intersect(InputData.triangles(mf(i),:),InputData.triangles(me(i),:));
        oftpa = setdiff(1:3,ia);
        oftpb = setdiff(1:3,ib);
        Bend = [Bend; [link,InputData.triangles(mf(i),oftpa),InputData.triangles(me(i),oftpb)]]; % [two nodes of bar, two nodes of remaining triangle points (from connected two InputData.triangles)]
    end
    PreprocessData.foldingHinges = Fold;
    PreprocessData.bendingHinges = Bend;
    
    pf0 = zeros(size(Fold,1),1); % Initial folding angle
    pb0 = zeros(size(Bend,1),1); % Initial bending angle
    for i = 1:size(Bend,1)
        pb0(i) = FoldKe(InputData.coordinates,Bend(i,:)); 
    end
    PreprocessData.initialFoldAngles = pf0;
    PreprocessData.initialBendAngles = pb0;
    
    % Generate bending stiffness array 
    bendingLengths = zeros(size(Bend,1),1); % Initialize bending hinge length array
    Kb = bendingLengths; % Initialize bending stiffness array
    for ang=1:size(Bend,1)
        node1=InputData.coordinates(Bend(ang,1),:);   
        node2=InputData.coordinates(Bend(ang,2),:);
        node3=InputData.coordinates(Bend(ang,3),:);   
        node4=InputData.coordinates(Bend(ang,4),:);
        bendingLengths(ang)=norm(node1-node2);
        A1=1/2*norm(cross(node2-node1,node3-node1));
        A2=1/2*norm(cross(node2-node1,node4-node1));
        Kb(ang)=0.1*1/6*bendingLengths(ang)^2*InputData.elasticModulus*InputData.thickness_sheet^3/(A1+A2);
    end
    PreprocessData.Kb = Kb;
    
    % Folding stiffness details
    PreprocessData.valleyFolds = [];
    PreprocessData.mountainFolds = [];
    PreprocessData.Kf = [];

end