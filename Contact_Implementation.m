function [contactStiffness, contactForces] = Contact_Implementation(slaveNodes, masterNodes, slaveElements, penaltyFactor)
% CONTACT_IMPLEMENTATION Computes contact stiffness and forces between slave and master surfaces
% 
% Inputs:
%   slaveNodes     - Coordinates of slave nodes
%   masterNodes    - Coordinates of master nodes
%   slaveElements  - Connectivity of slave elements
%   penaltyFactor  - Penalty factor for contact enforcement
%
% Outputs:
%   contactStiffness - Global contact stiffness matrix
%   contactForces    - Reaction forces on contact surfaces (cell array: {slaveForces, masterForces})

global Gap
Gap = [ ];

% Initialize contact stiffness matrix
totalDOFs = 2 * (size(slaveNodes, 1) + size(masterNodes, 1));
contactStiffness = zeros(totalDOFs);

% Get contact surfaces and their global DOF IDs
[slaveSurface, slaveDOF_IDs] = getContactSurface(1);
[masterSurface, masterDOF_IDs] = getContactSurface(2);

% Initialize contact force storage
numSlaveNodes = size(slaveSurface, 1);
numMasterNodes = size(masterSurface, 1);
slaveForces = zeros(numSlaveNodes, 4);  % [dofX, dofY, forceX, forceY]
masterForces = zeros(numMasterNodes, 4);
slavePressures = zeros(numSlaveNodes, 4);
masterPressures = zeros(numMasterNodes, 4);

numSlaveElements = numSlaveNodes - 1;
numMasterElements = numMasterNodes - 1;

% Main contact processing loop - slave nodes vs master elements
for slaveNodeIdx = 1:numSlaveNodes
    % Get DOF IDs for current slave node
    slaveDOF = slaveDOF_IDs(2*slaveNodeIdx-1:2*slaveNodeIdx);
    
    % Initialize force storage for this slave node
    slaveForces(slaveNodeIdx, :) = [slaveDOF, 0, 0];
    slavePressures(slaveNodeIdx, :) = [slaveDOF, 0, 0];
    
    % Check contact with all master elements
    for masterElemIdx = 1:numMasterElements
        % Get DOF IDs for master element nodes
        masterNode1DOF = masterDOF_IDs(2*masterElemIdx-1:2*masterElemIdx);
        masterNode2DOF = masterDOF_IDs(2*(masterElemIdx+1)-1:2*(masterElemIdx+1));
        allDOFs = [slaveDOF, masterNode1DOF, masterNode2DOF];
        
        % Get nodal coordinates
        slaveX = slaveSurface(slaveNodeIdx, 1);
        slaveY = slaveSurface(slaveNodeIdx, 2);
        masterX1 = masterSurface(masterElemIdx, 1);
        masterX2 = masterSurface(masterElemIdx+1, 1);
        masterY1 = masterSurface(masterElemIdx, 2);
        masterY2 = masterSurface(masterElemIdx+1, 2);
        elementCoords = [slaveX, masterX1, masterX2; 
                         slaveY, masterY1, masterY2];
        
        %% Fast contact check - skip if clearly not in contact
        if abs(slaveX - masterX1) > 0.04/600
            zeroForce = zeros(6, 1);
            masterForces(masterElemIdx, :) = [masterNode1DOF, masterForces(masterElemIdx, 3:4)] + [0, 0, zeroForce(3:4)'];
            masterForces(masterElemIdx+1, :) = [masterNode2DOF, masterForces(masterElemIdx+1, 3:4)] + [0, 0, zeroForce(5:6)'];
            masterPressures(masterElemIdx, :) = masterForces(masterElemIdx, :);
            masterPressures(masterElemIdx+1, :) = masterForces(masterElemIdx+1, :);
            continue;
        end
        
        %% Solve contact for this slave node - master element pair
        [contactForce, contactStiff, gap] = Contact_solver(penaltyFactor, elementCoords, slaveNodeIdx);
        
        % Skip if no contact
        if isempty(contactStiff)
            zeroForce = zeros(6, 1);
            masterForces(masterElemIdx, :) = [masterNode1DOF, masterForces(masterElemIdx, 3:4)] + [0, 0, zeroForce(3:4)'];
            masterForces(masterElemIdx+1, :) = [masterNode2DOF, masterForces(masterElemIdx+1, 3:4)] + [0, 0, zeroForce(5:6)'];
            masterPressures(masterElemIdx, :) = masterForces(masterElemIdx, :);
            masterPressures(masterElemIdx+1, :) = masterForces(masterElemIdx+1, :);
            continue;
        end
        
        % Store gap value
        Gap(slaveNodeIdx, 1) = gap;
        
        %% Assemble contact stiffness matrix
        for row = 1:size(contactStiff, 1)
            for col = 1:size(contactStiff, 2)
                rowDOF = allDOFs(row);
                colDOF = allDOFs(col);
                contactStiffness(rowDOF, colDOF) = contactStiffness(rowDOF, colDOF) + contactStiff(row, col);
            end
        end
        
        %% Compute and store contact forces
        % Calculate master element length for force distribution
        tangentVector = elementCoords(:, 3) - elementCoords(:, 2);
        elementLength = norm(tangentVector);
        
        % Update slave node forces
        slaveForces(slaveNodeIdx, :) = [slaveDOF, slaveForces(slaveNodeIdx, 3:4)] + [0, 0, contactForce(1:2)'] * elementLength;
        
        % Update master node forces
        masterForces(masterElemIdx, :) = [masterNode1DOF, masterForces(masterElemIdx, 3:4)] + [0, 0, contactForce(3:4)'] * elementLength;
        masterForces(masterElemIdx+1, :) = [masterNode2DOF, masterForces(masterElemIdx+1, 3:4)] + [0, 0, contactForce(5:6)'] * elementLength;
        
        % Update pressure values (forces per unit length)
        slavePressures(slaveNodeIdx, :) = [slaveDOF, slavePressures(slaveNodeIdx, 3:4)] + [0, 0, contactForce(1:2)'];
        masterPressures(masterElemIdx, :) = [masterForces(masterElemIdx, 1:2), masterForces(masterElemIdx, 3:4) * elementLength];
        masterPressures(masterElemIdx+1, :) = [masterForces(masterElemIdx+1, 1:2), masterForces(masterElemIdx+1, 3:4) / elementLength];
    end
end

% Package contact forces for output
contactForces = {slaveForces, masterForces};

end