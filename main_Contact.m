%% ================== Macro-meso-micro-scale Contact Solver: SCHC FRAMEWORK ================
%  ---------------------------------------
% This is an open multiscale contact solver based on Truncation method and
% Finite Element approach. The provided framework can solve the 2D rough
% contact problems considering surface topography and deterministic
% textures. Users can embed other contact models into the micro-scale 
% contact solver for replacement.
%  ---------------------------------------

%% This is an example of 'Multiscale Ironing Problem'

%% ---------------- License/Citation Notice ------------------
% I) For direct use of this code (without modifications):
%    If you use this code as-is to solve other problems and publish the 
%    results (in any form), you should contact the authors for permission.

% II) For modified versions of this code:
%    If you adapt and improve this code (method) in your work, you just need 
%    to make appropriate citation.
%  ---------------------------------------

%  ---------------------------------------
% Reference: 'Tao Yang, Hanmin Peng*, Guoqing Wang, Xiongfeng Tang, Yi Zheng.
% A unified scale-coupling homogenized contact model for elasto-plastic
% rough surfaces with textures. Tribology International (2025)'
%  ---------------------------------------

%%  ----------Software Provider------------
% The software (i.e., MATLAB) is provided by Yi Zheng
% and his affiliated institution ( College of Biophotonics,
% South China Normal University, Guangzhou, 510631, People’s Republic of China).
%  ---------------------------------------

% Please feel free to contact us with any questions! 
%  - yang.tao@nuaa.edu.cn


clear;
clc;
close all;

load('RAND.mat')

%% ================ Pretreatment ================
% Material parameters
global Parameter randomFAI Penalty
Parameter = struct;
Parameter.t = 1;   % thickness of 2D-problem
Parameter.E = 200e6;     
Parameter.miu = 0.3;
Parameter.E_rou = 10*Parameter.E;
Parameter.fai = 20;
% Fractal parameters
Parameter.Dim = 2.5;
Parameter.G = 1e-8;
Parameter.n_min = 5;
Parameter.n_max = 6;


% Control parameters of the iron movement
NUM_Y = 20;  Y_start = 1e-5; Y_end = 0e-6;    % Compress
NUM_X = 100;  X_start = -0.02; X_end = 0.02;  % Drag
DIS_Y = (Y_end-Y_start)/NUM_Y*ones(1,NUM_Y);
DIS_X = (X_end-X_start)/NUM_X*ones(1,NUM_X);
DIS_IND = [zeros(1,length(DIS_Y)), DIS_X; DIS_Y, zeros(1,length(DIS_X))]';
num_steps = length(DIS_IND);

% Solver parameters
global Penalty bt
Penalty = 5e11;
bt = 1e20;     % A penalty-like parameter using for boundary constraint
max_iter = 50;
tolerance = 1e-2; % Convergence tolerance
damping_factor = 1; % Damping factor in the Newton-Raphson-based iteration

% ================ Elastic matrix ================
global D 
D = Parameter.E/(1-Parameter.miu^2)*[ 1 Parameter.miu 0; Parameter.miu 1 0; 0 0 (1-Parameter.miu)/2 ]; % Elastic matrix (plane stress problem)

% =============== Macro-scale configuration ===============
global Element_s  Node_s  Node_m Element_m Node Element
fname1 = 'plan_down_800.txt';   
[Node_s, Element_s] = Readmesh( fname1 );
fname2 = 'iron_large.txt';
[Node_m, Element_m] = Readmesh( fname2 );
Node_s(:,1) = Node_s(:,1) - 1.1e-5;
Node_m(:,2) = Node_m(:,2) - min(Node_m(:,2)) + Y_start;
Node = [Node_s; Node_m ];
Element = [Element_s;Element_m+size(Node_s,1)];
Node_s_ini = Node_s; Node_m_ini = Node_m;  Node_ini = Node;
L_slave = 0.04/800;

% ================ Meso-scale texture ================
global z_tex K_tex
H_tex = 1*1e-6;
L_tex = 0.04/5;
z_tex = @(x) H_tex*(cos(2*pi/L_tex * x)+1);
E_tex = 0.1*Parameter.E;
K_tex = @(x) E_tex * t * L_slave / ( H_tex*(cos(2*pi/L_tex * x)+1) );    % K = E*A/L


% % Plot the mesh
% patch( 'vertices', Node_s, 'faces', Element_s, 'facecolor', [.9, .9, .9] ); hold on
% MeshPlot(Node_m, Element_m);
% axis equal; axis on;

% ================= Determine boundaries =================
[side_left, side_right, side_up, side_down] = find_side_rectangle(Element_s, Node_s); 
side_sphere = find_side_circular(Element_m, Node_m, -0.02 + 0.005, 0.05+Y_start, 0.05);
[~, ~, side_control, ~] = find_side_rectangle(Element_m, Node_m); 
C1 = sortrows([side_up',Node_s(side_up,1)],  2);
side_up = C1(:,1);
C2 = sortrows([side_sphere',Node_m(side_sphere,1)],  2);
side_sphere = C2(:,1);
side_sphere = flip(side_sphere);

% ================= Determine contact boundaries =================
global boundry_slave  boundry_master boundry_slave_glo boundry_master_glo side_control_global
boundry_slave = side_up;
boundry_master = side_sphere;
boundry_slave_glo = boundry_slave;
boundry_master_glo = boundry_master + size(Node_s(:,1),1);
side_control_global = side_control + size(Node_s(:,1),1);

% Discretization of Meso-scale texture
global Z_tex_ini  
Z_tex_ini = z_tex(Node_s(boundry_slave,1));   % Initial texture
% K_tex = E_tex * t * L_slave / ( H_tex*(cos(2*pi/L_tex * Node_s(boundry_slave,1))+1) );  
% NOTE: In this case, the texture stiffness is assumed to be infinite, so the stiffness will not be calculated

% % Plot boundaries
% scatter(Node_s(boundry_slave,1), Node_s(boundry_slave,2), "filled", 'r') % Slave contact boundary
% scatter(Node_m(boundry_master,1), Node_m(boundry_master,2), "filled", 'blue')  % Master contact boundary
% scatter(Node_s(side_down,1), Node_s(side_down,2), "filled", 'black')    % Dirichlet boundary
% scatter(Node_m(side_control,1), Node_m(side_control,2), "filled", 'black')    % Dirichlet boundary
% scatter(Node_t(:,1), Node_t(:,2), "filled", 'green')    % Texture


% ================= Initialization =================
global F_linear F_con Gap
F_linear = zeros(2*length(Node(:,1)),1);  F_reaction = zeros(2*length(Node(:,1)),1);
F_con = zeros(2*length(Node),1);
U = zeros(2*length(Node),1);
dis = zeros(2*length(Node),1);
Node_step = cell(1,num_steps);
Node_s_step = cell(1,num_steps);
Node_m_step = cell(1,num_steps);
dis_step = cell(1,num_steps);
dis_s_step = cell(1,num_steps);
dis_m_step = cell(1,num_steps);
F_con_slave = cell(1,num_steps);
GAP = cell(1,num_steps);
F_con_sum = zeros(length(num_steps),1);
[K_s_ini, ~, EA_s] = Linear_Elasticity(Node_s, Element_s, D);
[K_m_ini, ~, EA_m] = Linear_Elasticity(Node_m, Element_m, D);
[~, ~, EA] = Linear_Elasticity(Node, Element, D);
K_e_ini = [K_s_ini,zeros(size(K_s_ini, 1),size(K_m_ini,2)) ;zeros(size(K_m_ini, 1),size(K_s_ini,2)) ,K_m_ini];

%% ================= Solve the problem =================
tic   % Start timing
for step = 1:num_steps
    Node_m = Node_m + DIS_IND(step,:);
    Node_current = [Node_s; Node_m];
    Node = Node_current;
    residual = inf;   % Residual error
    node_delta = Node - Node_ini;
    U_iter = reshape(node_delta', 1, [])';
   
    fixed_dofs = [2*side_down-1, 2*side_down,2*side_control_global, 2*side_control_global-1];
    free_dofs = setdiff(1:length(U), fixed_dofs);

    %% ======== Newton-Raphson-based iteration ========
    iter = 0;
    while residual > tolerance && iter < max_iter
        % Initialize the stiffness matrix and residual vector of the current iteration (global)
        K_Global = zeros(2*length(Node), 2*length(Node));  % Global stiffness matrix
        R = zeros(2*length(Node), 1);        
        F_con = zeros(2*length(Node),1);
        F_con_s = zeros(2*length(Node),1);

        %% Contact part (penalty method)
        [K_Contact, RC_force] = Contact_Implementation(Node_s, Node_m, Element, Penalty);
        RC_force_slave = RC_force{1};     RC_force_master = RC_force{2};
        for i = 1:length(RC_force_slave(:,1))
            F_con(RC_force_slave(i,1)) = F_con(RC_force_slave(i,1)) + RC_force_slave(i,3);
            F_con(RC_force_slave(i,2)) = F_con(RC_force_slave(i,2)) + RC_force_slave(i,4);
        end
        for i = 1:length(RC_force_master(:,1))
            F_con(RC_force_master(i,1)) = F_con(RC_force_master(i,1)) + RC_force_master(i,3);
            F_con(RC_force_master(i,2)) = F_con(RC_force_master(i,2)) + RC_force_master(i,4);
        end
        for i = 1:length(RC_force_slave(:,1))
            F_con_s(RC_force_slave(i,1)) = F_con_s(RC_force_slave(i,1)) + RC_force_slave(i,3);
            F_con_s(RC_force_slave(i,2)) = F_con_s(RC_force_slave(i,2)) + RC_force_slave(i,4);
        end
        

        % Add contact stiffness (tangential stiffness matrix) to global matrix
        K_Global = K_e_ini + K_Contact;

        %% Boundary condition treatment
        K_Global(fixed_dofs, :) = 0;
        K_Global(:, fixed_dofs) = 0;
        K_Global(fixed_dofs, fixed_dofs) = eye(length(fixed_dofs)) * bt;
        

        % Residual force
        R_in = K_e_ini * U_iter;   % Equivalent to the internal force obtained by traversing each element
        R = R_in - F_con;
        
        R(fixed_dofs) = 0; % Residual zeroing of fixed degrees of freedom

        %% Solving linear equations: incremental displacement
        delta_U = zeros(2*length(Node), 1); 
        delta_U(free_dofs) = damping_factor * K_Global(free_dofs,free_dofs) \ (-R(free_dofs)); 

        %% Update displacement
        U_iter = U_iter + delta_U;
            Node = Node_ini + reshape(U_iter,2,length(U_iter)/2)';
            Node_s = Node(1:length(Node_s),:);
            Node_m = Node(length(Node_s)+1:end,:);
        
        %% Calculate residual L2 norm
        if norm(R_in) < 1e-3
            residual = norm(R);
        else
            residual = norm(R)/norm(R_in);
        end
        
        iter = iter + 1;
        RES(iter,step) = residual;
    end

    % Store the convergence displacement and other results
    U = U_iter;
    Node_step{step} = Node;
    Node_s_step{step} = Node_s;
    Node_m_step{step} = Node_m;
    dis_step{step} = U_iter;
    dis_s_step{step} = U_iter(1:2*length(Node_s));
    dis_m_step{step} = U_iter(2*length(Node_s)+1:end);
    F_con_slave{step} = RC_force_slave;
    F_con_sum(step) = abs(sum(F_con_s));
    GAP{step} = Gap;

    fprintf('Iteration %d: iter_num = %d  Res = %.4f\n', step-1, iter, residual);

% pause(0.2)
% figure
% MeshPlot(Node_step{step}, Element)

end
toc   % End timing


%% ================= Post-processing =================

%% Stress
sgm_step = cell(1,num_steps);
for i=1:num_steps
    for j=1:length(Element_s(:,1))
        sgm_s_step{i}(:,j) = stress_calculate(Element_s, Node_s_step{i}, j, EA_s, D, Parameter.t, dis_s_step{i});   
    end
    for j=1:length(Element_m(:,1))
        sgm_m_step{i}(:,j) = stress_calculate(Element_m, Node_m_step{i}, j, EA_m, D, Parameter.t, dis_m_step{i});  
    end
end


%% Transformation of nodal force into pressure  (p = Fc/An)
for i = 1:num_steps
    for j=1:length(Node_ini)
        node_glo(2*j-1) = Node_ini(j,1);
        node_glo(2*j) = Node_ini(j,2);
    end
    Cor_x = node_glo(F_con_slave{i}(:,1));   Cor_y = node_glo(F_con_slave{i}(:,2));
    force_x = F_con_slave{i}(:,3);  force_y = F_con_slave{i}(:,4);

    for j = 1:length(force_y)
        if j < length(force_y)
            A_seg = sqrt( (Cor_x(j+1)-Cor_x(j))^2  +  (Cor_y(j+1)-Cor_y(j))^2 ) * Parameter.t;
        else
            A_seg = sqrt( (Cor_x(j-1)-Cor_x(j))^2  +  (Cor_y(j-1)-Cor_y(j))^2 ) * Parameter.t;
        end
        pressure_x(j,i) = force_x(j)/A_seg;
        pressure_y(j,i) = force_y(j)/A_seg;
    end
end

% Plot the pressure curves of all loading steps
figure; hold on
for i=1:1:num_steps
    plot(Cor_x, -pressure_y(:,i), 'LineWidth', 2,'Color','blue');
    title('Contact pressure distribution');
    box on
end


%% Draw and save GIF
figure('Position', [150, 150, 1200, 900]*0.8); % [x, y, width, height]
filename = 'Body_stress.gif';
for i = 1:num_steps
    clf;
  
    % Body stress
%     hold on
%     Plot_stress(Element_s,Node_s_step{i}, sgm_s_step{i}, 4);     % Plot_stress(node, sgm, type)  4: von Mises
%     Plot_stress(Element_m,Node_m_step{i}, sgm_m_step{i}, 4);     % Plot_stress(node, sgm, type)  4: von Mises
%     axis([-0.03 0.03 -0.025 0.01]*1)
%     title('Von Mises Stress');
%     text(0,  0.005, sprintf('step = %g', i),'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'black');
%     axis off
%     caxis([0, 2.5]*1e7*1)

    % Deformed mesh
    hold on
    MeshPlot(Node_s_step{i}, Element_s)
    MeshPlot(Node_m_step{i}, Element_m)

    % Pressure distribution
%     hold on
%     plot(Cor_x, -pressure_y(:,i), 'LineWidth', 2,'Color','blue');
%     title('Contact pressure distribution');
%     box on

    % Surface displacement
%     hold on
%     disp = Node_s_step{i} - Node_s_ini;
%     plot(Node_s_ini(boundry_slave(:),1),disp(boundry_slave(:),1),'Color','red','LineWidth',2)
%     plot(Node_s_ini(boundry_slave(:),1),disp(boundry_slave(:),2),'LineStyle','--','Color','blue','LineWidth',2)
%     title('Surface displacement');
%     box on
%     legend('u_x', 'u_y', 'Location', 'best'); % 'Location', 'best' 

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end



% Save the results
DATA = struct;
DATA.Cor_x = Cor_x;
DATA.Element_m = Element_m;
DATA.Element_s = Element_s;
DATA.GAP = GAP;
DATA.Node_m_step = Node_m_step;
DATA.Node_s_step = Node_s_step;
DATA.Pressure_x = pressure_x;
DATA.Pressure_y = pressure_y;
DATA.sgm_m_step = sgm_m_step;
DATA.sgm_s_step = sgm_s_step;
DATA.fai = Parameter.fai;
DATA.h = H_tex;
DATA.RES = RES;



%%  ====================== Postprocessing ============================
% Note:Need to download a plugin called 'othercolor'

load('DATA.mat')   % Saved from the above calculation
Cor_x = DATA.Cor_x;
Element_m = DATA.Element_m;
Element_s = DATA.Element_s;
GAP = DATA.GAP;
Node_m_step = DATA.Node_m_step;
Node_s_step = DATA.Node_s_step;
pressure_y = DATA.Pressure_y;
sgm_m_step = DATA.sgm_m_step;
sgm_s_step = DATA.sgm_s_step;

num_steps = length(pressure_y(1,:));
load('Node_s_ini.mat')

%% Pressure distribution
plotPressureDistribution(Cor_x', -pressure_y');
axis([-0.02 0.02 0 120 0 4e5])

%% Normal contact force
L_slave = 0.04/800;
for i=1:num_steps
    F_sum(i) = -sum(pressure_y(:,i))*L_slave;
end
figure
plot(F_sum)

%% Normal displacement
figure; hold on
scale_factor = 2000;
for i=1:num_steps-12
    Node = Node_s_ini + (Node_s_step{i}(:,2) - Node_s_ini(:,2)) * scale_factor;
    Node(:,3) = 0.1 - i*0.1/num_steps;
    Node(:,3) = i*0.1/num_steps;
    NODE = [Node(:,3) Node(:,1) Node(:,2)];
    u = (Node_s_step{i} - Node_s_ini) * scale_factor;
    patch('Faces',Element_s,'Vertices',NODE,'facevertexcdata',u(:,2), 'facecolor','interp','EdgeColor', 'none'	)
    colormap(othercolor('RdYlBu11'));
    colormap(flipud(othercolor('RdYlBu11')));
    colorbar off
    axis equal
%     axis off
    set(gcf, 'Color', 'white');
    set(gca, 'FontSize', 14);
    caxis([-2 0.0]*1e-3)
    % caxis([-2 0.0]*1e-3 / scale_factor * 1e6)
    view([100, 300]);
    view(3);
end










