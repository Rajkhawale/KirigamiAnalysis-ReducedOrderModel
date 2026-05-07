function [PostprocessData, InputData] = Path_Analysis_BuildingFacade(InputData, PreprocessData)
tol = 3e-2; 
truss.U0 = zeros(3*size(InputData.coordinates,1),1);
U = truss.U0;

% Sheet loading increment size [rad]
f_inc = InputData.loadMagnitude/InputData.numberIncrements; 

% Fold and bending stiffness. Bar cross-sectional areas
angles.Kf = PreprocessData.Kf; 
angles.Kb = PreprocessData.Kb; 
truss.A = PreprocessData.barAreas;  

% Constitutive models
angles.CMbend = @(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,0,360);
angles.CMfold = @(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,0,360);
truss.CM = @(Ex)Ogden(Ex, InputData.elasticModulus); 

% Geometry (nodes, bar) information
angles.fold = PreprocessData.foldingHinges;
angles.bend = PreprocessData.bendingHinges;
truss.Bars = InputData.bar_elems;
truss.Node = InputData.coordinates;
truss.B = PreprocessData.B; 
truss.L = PreprocessData.L_bar;

% Initial angles
angles.pf0 = PreprocessData.initialFoldAngles;
angles.pb0 = PreprocessData.initialBendAngles;

truss.FixedDofs = PreprocessData.fixedDOFs; 


if strcmpi(InputData.loadType, 'Force')
    MaxIcr = InputData.numberIncrements;
    b_lambda = InputData.loadMagnitude/InputData.numberIncrements;
    % b_lambda = 1/InputData.numberIncrements;
    Uhis = zeros(3*size(InputData.coordinates,1),MaxIcr);
    FreeDofs = setdiff(1:3*size(InputData.coordinates,1),truss.FixedDofs);
    lmd = 0; icrm = 0; MUL = [U,U];
    lamHis = zeros(InputData.numberIncrements,1);
    angles.pf0_Orig=angles.pf0;
    F = PreprocessData.DOF_load;
    while icrm<MaxIcr && lmd < 1 % Analysis stops when applied load is below assigned load (\lambda < 1)
        icrm = icrm+1;
        %% Equivalent decrement/increment in angles
        if InputData.foldStep
            angles.pf0(PreprocessData.mountainFolds)=angles.pf0_Orig(PreprocessData.mountainFolds)-f_inc*icrm;
            angles.pf0(PreprocessData.valleyFolds)=angles.pf0_Orig(PreprocessData.valleyFolds)+f_inc*icrm;
        end
        iter = 0; err = 1;
        fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
        while err>tol && iter<InputData.maxIterations
            iter = iter+1; 
            [IF,K] = GlobalK_fast_ver(U,InputData.coordinates,truss,angles);
            R = lmd*F-IF;   
            MRS = [F,R];
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
            dUp = MUL(:,1); 
            dUr = MUL(:,2);
            if iter==1
                dUr = 0*dUr; 
            end
            % Modified Generalized Displacement Control Method]
            if iter==1
                if icrm==1
                    sinal=sign(dot(dUp,dUp));
                    dlmd=b_lambda;
                    numgsp=dot(dUp,dUp);   
                else
                    sinal=sinal*sign(dot(dupp1,dUp));
                    gsp=numgsp/dot(dUp,dUp);
                    dlmd=sinal*b_lambda*sqrt(gsp);
                end 
                dupp1=dUp;
                dupc1=dUp;
            else
                dlmd=-dot(dupc1,dUr)/dot(dupc1,dUp);
            end
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            err = norm(dUt(FreeDofs));
            lmd = lmd+dlmd;
            fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8
                disp('Divergence!')
                break
            end
        end

        if iter>15
            b_lambda = b_lambda/2;
            disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = lamHis(max(icrm,1));   % restore load
        elseif iter<3
            disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            lamHis(icrm) = lmd; 
        else
            Uhis(:,icrm) = U;
            lamHis(icrm) = lmd; 
        end
    end

elseif strcmpi(InputData.loadType, 'Displacement')
    Uhis = zeros(3*size(InputData.coordinates,1),InputData.numberIncrements*2);
    % A_CutOpen = zeros(21,2); % As we are storing total cut open area for each 1% stretching and our total stretching is 20%.
    % A_CutOpen(:,1) = [0:1/20:1];
    % flag = 2; % This flag for picking appropriate indices from A_CutOpen when the stretching complete 1%.
    Fdsp = PreprocessData.DOF_load/InputData.numberIncrements;
    ImpDofs = find(Fdsp~=0);
    FreeDofs = setdiff(setdiff(1:3*size(InputData.coordinates,1),truss.FixedDofs),ImpDofs);
    icrm = 0;  
    dspmvd = 0;  
    attmpts = 0;
    mvstepsize = 1;  
    damping = 1;
    % To get RF for fixed nodes
    ImpDofs_fixed = truss.FixedDofs;

    % To get degree of freedoms for reaction force plot nodes
    InputData.RF_dofs = 3*InputData.nodes_bottom-1; % For Y direction nodes
    Fhis = zeros(numel(InputData.RF_dofs), InputData.numberIncrements); 

    angles.pf0_Orig=angles.pf0;
    while dspmvd <= 1 && attmpts <= 20      
        icrm = icrm+1;
        %% Equivalent decrement/increment in angles
        if InputData.foldStep
            angles.pf0(PreprocessData.mountainFolds)=angles.pf0_Orig(PreprocessData.mountainFolds)-f_inc*icrm;
            angles.pf0(PreprocessData.valleyFolds)=angles.pf0_Orig(PreprocessData.valleyFolds)+f_inc*icrm;
        end
        iter = 0; err = 1;   
        fprintf('icrm = %d, dspimps = %6.4f\n',icrm,dspmvd);
        U = U+mvstepsize*Fdsp;
        U(truss.FixedDofs)=0;
        while err>tol && iter<InputData.maxIterations
            iter = iter+1;
            [IF,K] = GlobalK_fast_ver(U,InputData.coordinates,truss,angles);
            
            % For my mesh 
            % load_node = ceil(((InputData.mesh_scale_x+1)*(InputData.mesh_scale_y+1))/2);

            % For Matlab generated Mesh
            % target_point = [0.03,-0.07]; % For circle with two cut 
            target_point = [0,0.18]; % For circle with one cut
            dists = vecnorm(InputData.coordinates(:, 1:2) - target_point, 2, 2);
            [~, load_node] = min(dists);

            target_point2 = [0,0.55];% [0,-0.05]; % For square with one cut
            dists2 = vecnorm(InputData.coordinates(:, 1:2) - target_point2, 2, 2);
            [~, load_node2] = min(dists2);

            F = zeros(length(U),1);
            % F(load_node*3) = 1e-5; 
            % F(load_node2*3) = 1e-5;

            % if icrm< 400
            %     F(load_node*3) = 1e-50;  % Adding out-of-plane load at the center of plate. 
            % else
            %     F(load_node*3) = 0;
            % end

            IF = full(IF);
            IF = IF+F;
            IF = sparse(IF);
            dU = zeros(3*size(InputData.coordinates,1),1);
            dU(FreeDofs) = K(FreeDofs,FreeDofs)\(-IF(FreeDofs));
            err = norm(dU(FreeDofs));
            U = U+damping*dU; 
            fprintf('iter = %d, err = %6.4f, icrm = %d, attempts = %d, dspmvd = %6.4f\n',iter,err,icrm, attmpts, dspmvd);
        end

        if iter>=((mvstepsize>1)+1)*InputData.maxIterations/(damping+1)  
            % an aggressive step needs more iterations
            attmpts = attmpts+1;
            icrm = icrm-1;
            if attmpts<=20  
                mvstepsize = mvstepsize*0.5; 
                disp('Take a more conservative step...')
            % else
            %     mvstepsize = max(mvstepsize,1)*1.5;  
            %     damping = damping*0.75;
            %     disp('Take a more aggressive step...')
            end
            U = Uhis(:,max(icrm,1)); % restore displacement            
        else
            
            Uhis(:,icrm) = U;
            [Fend,~] = GlobalK_fast_ver(U,InputData.coordinates,truss,angles);
            Fhis(:,icrm) = Fend(InputData.RF_dofs)'; 

            % Computing and storing cut opening area information at specific interval. (after each 1% of stretching)
            % if  dspmvd > A_CutOpen(flag,1) || dspmvd > (1-1/InputData.numberIncrements) % The last statement consideres last stretching increment.
            %     Ui = reshape(Uhis(:,icrm), 3, [])';  % [X Y Z] displacements, size: num_nodes × 3 % Extract displacement at current increment
            %     deformed_coords = InputData.coordinates(:,1:3) + Ui; % New positions = original + displacement
            %     A_open = zeros(length(InputData.cut_node_ids),1);
            % 
            %     for i = 1:length(InputData.cut_node_ids)
            %         ids = InputData.cut_node_ids{i};
            %         coords_def = deformed_coords(ids, 1:2);  % Only (X, Y) for projected area
            %         % Compute polygon area from boundary
            %         k = convhull(coords_def(:,1), coords_def(:,2));
            %         A_open(i) = polyarea(coords_def(k,1), coords_def(k,2));
            %     end
            %     A_CutOpen(flag,2) = sum(A_open);
            %     flag = flag+1;
            %     fprintf('flag= %d\n',flag);
            % end

            % Resetting the updating parameters
            dspmvd = dspmvd+mvstepsize/InputData.numberIncrements;
            attmpts = 0; % Resetting the attmpts
            damping = 1;  % Resetting the damping
            if mvstepsize<1
                mvstepsize = min(mvstepsize*1.1,1); % gradually go back to 1
            else
                mvstepsize = max(mvstepsize*0.9,1);
            end
        end

        

    end
else
    disp('Unknown load type!!!')
end

if strcmpi(InputData.loadType,'Force') % lamHis gives the proportion of the load, not the actual load
    Fhis = lamHis*nonzeros(F)';
    lamHis(icrm+1:end,:) = [];
    PostprocessData.forceAchieved = [InputData.loads(:,1), lamHis(end)*InputData.loads(:,2:4)];
end

Uhis(:,icrm+1:end) = [];
% Fhis(icrm+1:end,:) = [];

PostprocessData.Uhis = real(Uhis);
PostprocessData.Fhis = real(Fhis);
% PostprocessData.A_CutOpen = A_CutOpen;
PostprocessData.foldAngles = angles.pf0;
if strcmpi(InputData.loadType,'Force')
    PostprocessData.lamHis = lamHis;
end
end
