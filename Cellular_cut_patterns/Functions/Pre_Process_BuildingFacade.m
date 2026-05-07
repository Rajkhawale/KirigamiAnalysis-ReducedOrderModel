

function [PreprocessData, InputData] = Pre_Process_BuildingFacade(InputData)

    %% Z Perturbation
    % 
    % max_disp = 0.05;             % Maximum Z displacement
    % x_center = 1;               % X center of cut
    % y_cut = 1.2133;               % Y location of the top of the cut
    % y_max = max(InputData.coordinates(:,2));  % Top of circular sheet
    % 
    % % Get all nodes above the cut
    % coords = InputData.coordinates(:, 1:2);
    % nodes_above_cut = find(coords(:,2) >= y_cut);
    % 
    % % Extract X and Y of those nodes
    % x_vals = coords(nodes_above_cut, 1);
    % y_vals = coords(nodes_above_cut, 2);
    % 
    % % Normalize distances from cut center (for Gaussian-like decay)
    % x_spread = 1.0;  % Controls lateral decay width
    % y_spread = y_max - y_cut;  % Full vertical span
    % 
    % % Distance from centerline (in x)
    % dx = (x_vals - x_center) / x_spread;
    % dy = (y_vals - y_cut) / y_spread;
    % 
    % % 2D Gaussian weights: centered at (x=0, y=y_cut)
    % weights = exp(-dx.^2 - dy.^2);
    % 
    % % Apply Z perturbation
    % InputData.coordinates(nodes_above_cut, 3) = max_disp * weights;

    % Maximum Z displacement at center
    max_disp = InputData.perturb;

    % Center of the cut (x and y)
    x_center = 1;
    y_center = 1.2;

    % Coordinates
    coords = InputData.coordinates(:, 1:2);
    x_vals = coords(:, 1);
    y_vals = coords(:, 2);

    % Spreads (controls how wide the Gaussian decay is)
    x_spread = 1.0;
    y_spread = 1.0;  % Use a consistent value, or adjust based on your geometry

    % Normalized distances from center
    dx = (x_vals - x_center) / x_spread;
    dy = (y_vals - y_center) / y_spread;

    % 2D Gaussian weights: centered at (x_center, y_center)
    weights = exp(-dx.^2 - dy.^2);

    % Apply Z perturbation
    InputData.coordinates(:, 3) = max_disp * weights;

    %% Works well for 5 horz and 10 vert cells scenario. The curection of opening of the cuts get correct with this perturbation.
    % Maximum Z displacement
    % max_disp = -0.02;
    % 
    % % Coordinates
    % coords = InputData.coordinates(:, 1:2);
    % y_vals = coords(:, 2);  % Extract Y coordinates
    % 
    % % Y range
    % y_min = 0;
    % y_max = 2.4;
    % 
    % % Normalize Y distance (0 at bottom → 1 at top)
    % dy = (y_vals - y_min) / (y_max - y_min);
    % 
    % % 1D Gaussian-like decay from bottom (y = 0) to top (y = 2.4)
    % weights = exp(-dy.^2);  % max at bottom (dy=0), decay toward top
    % 
    % % Apply Z perturbation
    % InputData.coordinates(:, 3) = max_disp * weights;



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
    
    tol = 0.095; % Define a small tolerance to identify edge-aligned nodes
    y_coords = InputData.coordinates(:,2);
    y_min = min(y_coords); % Get Y limits of the rectangle
    y_max = max(y_coords);
    InputData.nodes_bottom = find(abs(y_coords - y_min) < tol);
    InputData.nodes_top = find(abs(y_coords - y_max) < tol);

    InputData.supports = [InputData.nodes_bottom, ones(length(InputData.nodes_bottom),1), ones(length(InputData.nodes_bottom),1), ones(length(InputData.nodes_bottom),1); % Fixing in all direction
                          InputData.nodes_top, ones(length(InputData.nodes_top),1), zeros(length(InputData.nodes_top),1), ones(length(InputData.nodes_top),1)]; % Fixing in X & Z direction
    
    InputData.loads = [InputData.nodes_top, zeros(length(InputData.nodes_top),1), InputData.loadMagnitude*ones(length(InputData.nodes_top),1), zeros(length(InputData.nodes_top),1)]; % Stretching Y upward
    
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
    %Shear coefficients
% a0 = 0.25; a1 = 0.1; a2 = 0.25; % Original values -> a0 = 0.45; a1 = 0; a2 = 0.4;
a0 = 0.25; a1 = 0; a2 = 0.5; 

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
        % if L(j) == Lengths(3)
        %     A(j) = (a2*AspectRatio(i)^2 + a1*AspectRatio(i) + a0)*panelArea*InputData.thickness_sheet/(L(j)); 
        % else
        %     % A(j) = 0.3*panelArea*InputData.thickness_sheet/(L(j)); % original factor = 0.5
        %     A(j) = (a2*AspectRatio(i)^2 + a1*AspectRatio(i) + a0)*panelArea*InputData.thickness_sheet/(L(j)); 
        % end
        if L(j) == Lengths(3)
            A(j) = panelArea*InputData.thickness_sheet/(L(j)); 
        else
            % A(j) = 0.3*panelArea*InputData.thickness_sheet/(L(j)); % original factor = 0.5
            A(j) = panelArea*InputData.thickness_sheet/(L(j)); 
        end
    end; clear j

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
        Kb(ang)=1*1/6*bendingLengths(ang)^2*InputData.elasticModulus*InputData.thickness_sheet^3/(A1+A2); % original factor = 1
    end
    PreprocessData.Kb = Kb;
    
    % Folding stiffness details
    PreprocessData.valleyFolds = [];
    PreprocessData.mountainFolds = [];
    PreprocessData.Kf = [];

end

