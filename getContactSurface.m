function [P, Global_DOF_ID] = getContactSurface(Slave_or_Master)
% getContactSurface retrieves the contact surface coordinates and the corresponding
% global degrees of freedom (DOF) indices for either the slave or master surface.

% Input:
%   Slave_or_Master - 1 for Slave surface, 2 for Master surface

% Output:
%   P - The coordinates of the contact surface nodes
%   Global_DOF_ID - The global degrees of freedom (DOF) indices for the contact surface

global boundry_master_glo boundry_slave_glo Node

switch Slave_or_Master
    case 2  % Master surface
        P = Node(boundry_master_glo, :); % Coordinates of the master contact surface nodes
        Global_Node_ID = boundry_master_glo;

        % Initialize the global DOF indices for the master surface
        for i = 1:length(Global_Node_ID)
            Global_DOF_ID(2*i-1) = 2 * Global_Node_ID(i) - 1; % Global DOF index for x-coordinate
            Global_DOF_ID(2*i) = 2 * Global_Node_ID(i);         % Global DOF index for y-coordinate
        end

    case 1  % Slave surface
        P = Node(boundry_slave_glo, :); % Coordinates of the slave contact surface nodes
        Global_Node_ID = boundry_slave_glo;

        % Initialize the global DOF indices for the slave surface
        for i = 1:length(Global_Node_ID)
            Global_DOF_ID(2*i-1) = 2 * Global_Node_ID(i) - 1; % Global DOF index for x-coordinate
            Global_DOF_ID(2*i) = 2 * Global_Node_ID(i);         % Global DOF index for y-coordinate
        end

end

end
