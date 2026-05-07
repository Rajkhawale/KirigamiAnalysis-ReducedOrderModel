function PostprocessData = Post_Process_Kirigami(InputData, PreprocessData, PostprocessData)

barConstitutiveModel = @(Ex)Ogden(Ex, InputData.elasticModulus);
bendConstitutiveModel = @(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,0,360);
foldConstitutiveModel = @(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,0,360);

%% Get Data
Exbar = zeros(size(InputData.bar_elems,1),size(PostprocessData.Uhis,2)); 
FdAngle = zeros(size(PreprocessData.foldingHinges,1),size(PostprocessData.Uhis,2)); 
BdAngle = zeros(size(PreprocessData.bendingHinges,1),size(PostprocessData.Uhis,2));
for icrm=1:size(PostprocessData.Uhis,2)
    Ui = PostprocessData.Uhis(:,icrm);
    Nodenw = zeros(size(InputData.coordinates));
    Nodenw(:,1) = InputData.coordinates(:,1)+Ui(1:3:end);
    Nodenw(:,2) = InputData.coordinates(:,2)+Ui(2:3:end);
    Nodenw(:,3) = InputData.coordinates(:,3)+Ui(3:3:end);
    
    eDofb = kron(InputData.bar_elems,3*ones(1,3))+repmat([-2,-1,0],size(InputData.bar_elems,1),2);
    du = Ui(eDofb(:,1:3))-Ui(eDofb(:,4:6));
    Exbar(:,icrm) = PreprocessData.B*Ui./PreprocessData.L_bar+0.5*sum(du.^2,2)./(PreprocessData.L_bar.^2);

    for del = 1:size(PreprocessData.bendingHinges,1)
        bend = PreprocessData.bendingHinges(del,:);
        BdAngle(del,icrm) = FoldKe(Nodenw,bend);
    end

    for fel = 1:size(PreprocessData.foldingHinges,1)
        fold = PreprocessData.foldingHinges(fel,:);
        FdAngle(fel,icrm) = FoldKe(Nodenw,fold);
    end
end

%% Interpret Data
[Sx_bar, ~, Wb] = barConstitutiveModel(Exbar);
Rspr_fd = zeros(size(FdAngle)); 
Efold = Rspr_fd;
Rspr_bd = zeros(size(BdAngle)); 
Ebend = Rspr_bd;
FdAngle_2 = FdAngle - pi*ones(size(FdAngle));
for i = 1:size(PostprocessData.Uhis,2)
    [Rspr_fdi, ~, Efoldi] = foldConstitutiveModel(abs(FdAngle_2(:,i)),...
                                                  abs(InputData.restAngleRad)*ones(size(PreprocessData.initialFoldAngles)),...
                                                  PreprocessData.Kf,...
                                                  PreprocessData.L_bar((size(PreprocessData.bendingHinges,1)+1):(size(PreprocessData.bendingHinges,1)+size(PreprocessData.foldingHinges,1))));
    [Rspr_bdi, ~, Ebendi] = bendConstitutiveModel(BdAngle(:,i),...
                                                  PreprocessData.initialBendAngles,...
                                                  PreprocessData.Kb,...
                                                  PreprocessData.L_bar(1:size(PreprocessData.bendingHinges,1)));
    Rspr_fd(:,i) = Rspr_fdi;
    Efold(:,i) = Efoldi;
    Rspr_bd(:,i) = Rspr_bdi; 
    Ebend(:,i) = Ebendi;
end; clear i

for i = 1:size(PostprocessData.Uhis,2)
    deformedDOFs = PostprocessData.Uhis(:,i);
                    deformedNodes = InputData.coordinates;
                    deformedNodes(:,1) = InputData.coordinates(:,1)+deformedDOFs(1:3:end);
                    deformedNodes(:,2) = InputData.coordinates(:,2)+deformedDOFs(2:3:end);
                    deformedNodes(:,3) = InputData.coordinates(:,3)+deformedDOFs(3:3:end); 
    PostprocessData.deformedNodes{i} = deformedNodes;
end; clear i

%% Collect results
    % Bars
    PostprocessData.barEnergy = diag(PreprocessData.L_bar.*PreprocessData.barAreas)*Wb;
    PostprocessData.barStrain = Exbar;
    PostprocessData.barStress = Sx_bar;
    % Folding hinges
    PostprocessData.foldingEnergy = Efold;
    PostprocessData.foldAngles = FdAngle;
    PostprocessData.foldResistantMoment = Rspr_fd;
    % Bending hinges
    PostprocessData.bendingAngle = BdAngle;
    PostprocessData.bendingResistantMoment = Rspr_bd;
    PostprocessData.bendingEnergy = Ebend;

end