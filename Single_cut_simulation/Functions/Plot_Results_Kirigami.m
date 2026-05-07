function Plot_Results_Kirigami(InputData, PreprocessData, PostprocessData)

    if ischar(InputData.plotIncrement)
        InputData.plotIncrement = size(PostprocessData.Uhis,2);
    end

    plotBarEnergy = PostprocessData.barEnergy(:,InputData.plotIncrement);
    plotFoldingEnergy = PostprocessData.foldingEnergy(:,InputData.plotIncrement);
    plotBendingEnergy = PostprocessData.bendingEnergy(:,InputData.plotIncrement);
    
    %% Sum total energy

    totalBarEnergy = sum(plotBarEnergy);
    totalFoldingEnergy = sum(plotFoldingEnergy);
    totalBendingEnergy = sum(plotBendingEnergy);

    %% Deformed and Bar Energy Plot

    if strcmpi(InputData.plotDeformedShape,'yes')
            figure()
                deformedDOFs = PostprocessData.Uhis(:,InputData.plotIncrement);
                deformedNodes = InputData.coordinates;
                deformedNodes(:,1) = InputData.coordinates(:,1)+deformedDOFs(1:3:end);
                deformedNodes(:,2) = InputData.coordinates(:,2)+deformedDOFs(2:3:end);
                deformedNodes(:,3) = InputData.coordinates(:,3)+deformedDOFs(3:3:end); 

            % plot initial configuration
            if strcmpi(InputData.plotReferenceShape,'yes') 
                Bar_els = InputData.bar_elems;
                for link=1:size(Bar_els,1)
                    pts=InputData.coordinates(Bar_els(link,:),:);
                    plot3(pts(:,1),pts(:,2),pts(:,3),'Color',[165 164 164]./255,'Linewidth',1); hold on;
                end  
            end
            % plot panels 
            cellArray_triangles = mat2cell(InputData.triangles, ones(1, size(InputData.triangles, 1)), size(InputData.triangles, 2));
            for el=1:size(InputData.triangles,1)
                pts=deformedNodes(cellArray_triangles{el},:);
                patch(pts(:,1),pts(:,2),pts(:,3),[220, 220, 220]./255,'EdgeColor','none'); hold on;
            end

            % Plot truss bars
            Edge_els=setdiff(InputData.bar_elems,PreprocessData.bendingHinges(:,[1,2]),'rows');
            Bar_els=InputData.bar_elems; 
            Bar_els=setdiff(Bar_els,Edge_els,'rows');
                for link=1:size(Bar_els,1)
                    pts=deformedNodes(Bar_els(link,:),:);
                    % plot3(pts(:,1),pts(:,2),pts(:,3),'Color',[180, 180, 180]./255,'Linewidth',7.0); hold on;
                    plot3(pts(:,1),pts(:,2),pts(:,3),'Color',[0 0 0]./255,'Linewidth',1.3); hold on;
                end  

            % Plot Perimeter
            for link=1:size(Edge_els,1)
                pts=deformedNodes(Edge_els(link,:),:); % 81
                plot3(pts(:,1),pts(:,2),pts(:,3),'Color',[0 0 0]./255,'Linewidth',0.8); hold on;
            end    

            axis equal; axis off;
            % if strcmp(InputData.testType,'SineWave')
            %     view(57,20);
            % elseif strcmp(InputData.testType,'CurvedSquare')
            %     view(60,35);
            % elseif strcmp(InputData.testType,'FullAnnulus')
            %     view(15,21);
            % elseif strcmp(InputData.testType,'Canopy')
            %     view(40,25);
            % elseif strcmp(InputData.testType,'CutAnnulus')
            %     view(-15,45);
            % end
    end


    if strcmpi(InputData.plotEnergy,'yes')
        figure('Position',[100 100 1200 400])        
        subplot(1,3,1)
                patch('faces', InputData.triangles, 'vertices', InputData.coordinates, 'facecolor', [1 1 1]*0.95, ...
                      'linestyle', 'none', 'facelighting', 'flat', 'edgecolor', (0)*[1 1 1]); hold on;
                max_Energy=max(plotFoldingEnergy);
                min_Energy=min(plotFoldingEnergy);
                clr_range=max_Energy-min_Energy;
                clr_map=colormap('viridis');
                clr_num=size(clr_map,1);

                for el=1:size(plotFoldingEnergy,1)
                    nds=PreprocessData.foldingHinges(el,[1,2]);
                    if clr_range == 0
                        color_1 = 1;
                    else
                        color_1=ceil(clr_num*(plotFoldingEnergy(el)-min_Energy)/clr_range);
                        if color_1 == 0
                            color_1 = 1;
                        end
                    end
                    plot3(InputData.coordinates(nds,1), InputData.coordinates(nds,2), InputData.coordinates(nds,3),'Color',clr_map(color_1,:),'Linewidth',2)
                end

                colorbar('southoutside')
                caxis([0 max_Energy]);

                axis equal; axis off;
                light
                view(0,90)
                rotate3d on
                title('Fold Energy')
                set(gcf,'renderer','painters')

        subplot(1,3,2)
                patch('faces', InputData.triangles, 'vertices', InputData.coordinates, 'facecolor', [1 1 1]*0.95, ...
                      'linestyle', 'none', 'facelighting', 'flat', 'edgecolor', (0)*[1 1 1]); hold on;
                max_Energy=max(plotBarEnergy);
                min_Energy=min(plotBarEnergy);
                clr_range=max_Energy-min_Energy;
                clr_map=colormap('viridis');
                clr_num=size(clr_map,1);

                for el=1:size(plotBarEnergy,1)
                    nds=InputData.bar_elems(el,[1,2]);
                    if clr_range == 0
                        color_1 = 1;
                    else
                        color_1=ceil(clr_num*(plotBarEnergy(el)-min_Energy)/clr_range);
                        if color_1 == 0
                            color_1 = 1;
                        end
                    end
                    plot3(InputData.coordinates(nds,1), InputData.coordinates(nds,2),InputData.coordinates(nds,3),'Color',clr_map(color_1,:),'Linewidth',2)
                end

                colorbar('southoutside')
                caxis([0 max_Energy]);

                axis equal; axis off;
                light
                view(0,90)
                rotate3d on
                title('Bar Energy')
                set(gcf,'renderer','painters')
           subplot(1,3,3)
                patch('faces', InputData.triangles, 'vertices', InputData.coordinates, 'facecolor', [1 1 1]*0.95, ...
                'linestyle', 'none', 'facelighting', 'flat', 'edgecolor', (0)*[1 1 1]); hold on;
                max_Energy=max(plotBendingEnergy);
                min_Energy=min(plotBendingEnergy);
                clr_range=max_Energy-min_Energy;
                clr_map=colormap('viridis');
                clr_num=size(clr_map,1);

                for el=1:size(PreprocessData.bendingHinges,1)
                   nds=PreprocessData.bendingHinges(el,[1,2]);
                   if clr_range == 0
                       color_1 = 1;
                   else
                       color_1=ceil(clr_num*(plotBendingEnergy(el)-min_Energy)/clr_range);
                        if color_1 == 0
                            color_1 = 1;
                        end
                   end
                   plot3(InputData.coordinates(nds,1), InputData.coordinates(nds,2), InputData.coordinates(nds,3),'Color',clr_map(color_1,:),'Linewidth',2)
                end

                colorbar('southoutside')
                caxis([0 max_Energy]);

                axis equal; axis off;
                light
                view(0,90)
                rotate3d on
                title('Bend Energy')
    end           
end