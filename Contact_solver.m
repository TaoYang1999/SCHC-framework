function [contactForce, contactStiffness, normalGap] = Contact_solver(penalty, elementCoordinates, nodeIndex)
% ContactSolver calculates the contact force, contact stiffness, and normal gap for the given contact surface.
% 
% Input:
%   penalty             - Penalty factor for the contact interaction.
%   elementCoordinates  - Coordinates of the element nodes, containing both endpoints of the segment.
%   nodeIndex           - Index of the current node for texture gap modification.
%
% Output:
%   contactForce        - Contact force vector.
%   contactStiffness    - Contact stiffness matrix.
%   normalGap           - The normal gap function after texture correction.

global  Z_tex_ini   

% Constants for clarity
ZERO = 0.0; 
ONE = 1.0; 
EPSILON = 1e-6; 
SMALL_OFFSET = 0.00; 

% Initialize outputs
contactForce = []; 
contactStiffness = [];

% Calculate the tangential vector
tangentialVector = elementCoordinates(:,3) - elementCoordinates(:,2);
segmentLength = norm(tangentialVector);

% Check if the segment length is non-zero
if segmentLength < EPSILON
    return;
end

% Normalize the tangential vector and compute the normal vector
tangentialVector = tangentialVector / segmentLength;
normalVector = [-tangentialVector(2); tangentialVector(1)];

% Compute the normal gap function: Gn = (X_s - X_1) . N
normalGap = (elementCoordinates(:,1) - elementCoordinates(:,2))' * normalVector;

% Adjust normal gap based on texture information
textureGapCorrection = Z_tex_ini(nodeIndex);
normalGap = normalGap - textureGapCorrection;

% Check if the gap is positive (no contact)
% Uncomment the next line if needed
% if (normalGap >= ZERO), return; end

% Compute the natural coordinate at the contact point
alpha = (elementCoordinates(:,1) - elementCoordinates(:,2))' * tangentialVector / segmentLength;

% Check if the contact point is within the segment bounds
if (alpha > ONE + SMALL_OFFSET) || (alpha < -SMALL_OFFSET)
    return;
end

% Calculate the pressure at the contact point based on the gap
pressure = Pressure_cal_truncation(normalGap);

% If pressure is zero, return
if pressure == 0
    return;
end

% Debugging code, plot the contact force vector and contact pair (optional, can be removed later)
time = 4; % Set a fixed time for debugging
if time < 2
    scatter(elementCoordinates(1,1), elementCoordinates(2,1), 'filled', 'r');
    scatter(elementCoordinates(1,2:3), elementCoordinates(2,2:3), 'filled', 'g');
    theta = atan(normalVector(2) / normalVector(1));
    plot(elementCoordinates(1,2:3), elementCoordinates(2,2:3), 'b', 'LineWidth', 2);
    quiver(elementCoordinates(1,1), elementCoordinates(2,1), abs(normalGap) * cos(theta), abs(normalGap) * sin(theta), 'color', 'm', 'LineWidth', 2);
end

% Contact occurs in this segment, calculate the contact force multiplier
% lambda = -penalty * normalGap;   % Smooth surface
lambda = pressure; 

% Define the contact force direction vectors
normalDirection = [normalVector; -(ONE - alpha) * normalVector; -alpha * normalVector];
tangentialDirection = [ZERO; ZERO; -tangentialVector; tangentialVector];
contactDirection = normalDirection;

% Calculate the contact force
contactForce = lambda * contactDirection;

% Calculate the contact stiffness (unit stiffness on a singal contact pair)
contactStiffness = segmentLength * penalty * (contactDirection * contactDirection');
% Alternatively, the stiffness calculation could use the following:
% contactStiffness = segmentLength * (1 * pressure / abs(normalGap)) * (contactDirection * contactDirection');

end
