function [cleanFrame] = VoltImg_mapping_removeArtifact(currFrame)

% Parameters for setting FOV area that is gated
gateOff = 100; %pixel number of FOV in direction of scan gating goes off
gateOn = 400; %pixel number in FOV direction of scan where gating comes back on
    
% set the area of gated FOV area where statistics will be calculated
statArea = currFrame(:, gateOff:gateOn);
    
lineVar = zeros(size(statArea, 1), 1);
lineStdDev = zeros(size(statArea, 1), 1);
cleanFrame = currFrame;
    for frameLine = 1:size(statArea, 1)
        % Calculate statistics
        lineVar(frameLine) = var(statArea(frameLine, :));
        lineStdDev(frameLine) = std(statArea(frameLine, :));
        % NaN the lines that statistically indicates imaging artifact
        if lineVar(frameLine) > 5500
            cleanFrame(frameLine, :) = NaN;
        end
    end

end
     
        
        
