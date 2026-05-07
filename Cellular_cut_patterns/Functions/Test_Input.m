function [InputData, inputTest] = Test_Input(InputData)

    inputTest = 1;
           
    if ~isnumeric(InputData.restAngle)
        fprintf('ERROR: <restAngle> must be a number\n\n')
        inputTest = 0;
        restAngleRad = nan;
    elseif ~isreal(InputData.restAngle)
        fprintf('ERROR: <restAngle> must be a real number\n\n')
        inputTest = 0;
        restAngleRad = nan;
    elseif InputData.restAngle < -180 || InputData.restAngle > 180
            fprintf('ERROR: <restAngle> must be between -180 and 180 [deg]\n\n')
            inputTest = 0;
            restAngleRad = nan;
    elseif InputData.restAngle < 0
            restAngleRad = -(180 - abs(InputData.restAngle))*pi/180;
    else
            restAngleRad = (180 - InputData.restAngle)*pi/180;
    end
    
    if ~isnumeric(InputData.elasticModulus)
        fprintf('ERROR: <elasticModulus> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.elasticModulus)
        fprintf('ERROR: <elasticModulus> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.elasticModulus < 0
        fprintf('ERROR: <elasticModulus> must be greater than zero\n\n')
        inputTest = 0;
    end

    if ~isnumeric(InputData.thickness)
        fprintf('ERROR: <thickness> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.thickness)
        fprintf('ERROR: <thickness> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.thickness < 0
        fprintf('ERROR: <thickness> must be greater than zero\n\n')
        inputTest = 0;
    end

    if ~isnumeric(InputData.lengthScale)
        fprintf('ERROR: <lengthScale> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.lengthScale) 
        fprintf('ERROR: <lengthScale> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.lengthScale < 0
        fprintf('ERROR: <lengthScale> must be greater than zero\n\n')
        inputTest = 0;
    end
    
    if ~isnumeric(InputData.loadMagnitude)
        fprintf('ERROR: <loadMagnitude> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.loadMagnitude)
        fprintf('ERROR: <loadMagnitude> must be a real number\n\n')
        inputTest = 0;
    end

%     if ~ischar(InputData.testType)
%         fprintf('ERROR: <testType> must be one of the following:\n\n\t\tFullAnnulus\n\t\tSineWave\n\t\tCurvedSquare\n\t\tCanopy\n\t\tCutAnnulus\n\n')
%         inputTest = 0;
%     elseif ~strcmp(InputData.testType,'FullAnnulus') && ~strcmp(InputData.testType,'SineWave') && ~strcmp(InputData.testType,'CurvedSquare') && ~strcmp(InputData.testType,'Canopy') && ~strcmp(InputData.testType,'CutAnnulus')
%         fprintf('ERROR: <testType> must be one of the following:\n\n\t\tFullAnnulus\n\t\tSineWave\n\t\tCurvedSquare\n\t\tCanopy\n\t\tCutAnnulus\n\n')
%         inputTest = 0;
%     end

    if ~isnumeric(InputData.numberDivisions) || ceil(InputData.numberDivisions) ~= InputData.numberDivisions
        fprintf('ERROR: <numberDivisions> must be an integer\n\n')
        inputTest = 0;
    elseif InputData.numberDivisions < 0
        fprintf('ERROR: <numberDivisions> must be a positive integer\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotNodes)
        fprintf('ERROR: <plotNodes> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotNodes,'yes') && ~strcmpi(InputData.plotNodes,'no')
        fprintf('ERROR: <plotNodes> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end
 
    if ~ischar(InputData.plotEnergy)
        fprintf('ERROR: <plotEnergy> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotEnergy,'yes') && ~strcmpi(InputData.plotEnergy,'no')
        fprintf('ERROR: <plotEnergy> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotDeformedShape)
        fprintf('ERROR: <plotDeformedShape> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotDeformedShape,'yes') && ~strcmpi(InputData.plotDeformedShape,'no')
        fprintf('ERROR: <plotDeformedShape> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotReferenceShape)
        fprintf('ERROR: <plotReferenceShape> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotReferenceShape,'yes') && ~strcmpi(InputData.plotReferenceShape,'no')
        fprintf('ERROR: <plotReferenceShape> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end
    
    if ischar(InputData.plotIncrement)
        if ~strcmpi(InputData.plotIncrement,'end')
            fprintf('ERROR: <plotIncrement> must be an integer between 1 and the specified <incrementNumber> or "end"\n\n')
            inputTest = 0;
        end
    elseif InputData.plotIncrement > InputData.numberIncrements
        fprintf('ERROR: <plotIncrement> might exceed available number of increments; proceed with caution or choose "end" instead\n\n')
        if inputTest == 1
            inputTest = 2;
        end
    end
            
    InputData.restAngleRad = restAngleRad;
 
    
           
    
    
    
