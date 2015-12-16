% Geometric parameters
geometry = struct();
geometry.DCB_height = 5e-3;
geometry.DCB_length = 100e-3;
geometry.DCB_width = 25e-3;
geometry.DCB_a0 = 25e-3;
geometry.DCB_elem_size_height = geometry.DCB_height/8; % 8 elements in the height
geometry.DCB_elem_size_length = geometry.DCB_a0/25; % 15 elements in the crack

% Load parameters
load_d = 5e-3; % Prescribed displacement

% Material paramters
material = struct();
material.E = 40e9;
material.poisson = 0.29;
material.penalty_stiffness = 1e14;
material.Gc = 613; % Critical energy release rate
material.tau0 = 12e6; % Onset traction

% Define materials
materialOptions = struct();
materialOptions(1).matID = 1;
materialOptions(1).options = [1, material.E, material.poisson]; % Linear elastic, Isotropic, Youngs modulus, Poisson ration
materialOptions(2).matID = 1;
materialOptions(2).options = [material.penalty_stiffness, material.Gc, material.tau0]; % This is for the cohesive zone

% Define elements
elementOptions = struct();
elementOptions(1).elemID = 1;
elementOptions(1).options = [geometry.DCB_width, 2, 1];  % Q4isoparametric, Thickness, Integration order, Plane stress
elementOptions(2).elemID = 3; % 2 == linear cohesive elements, 3 == mpc approach
elementOptions(2).options = [geometry.DCB_width, 2]; % Q4cohesive, Thickness, Integration order

% Solver options
solverOptions = struct();
solverOptions.perb = 0; % Pertubation for perturbed Lagrangian
solverOptions.tol = 1e-8; % This is compared with the norm of the residual
solverOptions.max_iter = 5e3; % This is the maximum allowed cumulative iterations
solverOptions.stepsize = 1/100; % Initial stepsize
solverOptions.min_stepsize = 1/1e4; % Minimum stepsize
solverOptions.max_stepsize = 1/100; % Maximum stepsize

% If using the MPC approach the mesh should add the MPCs 
if elementOptions(2).elemID == 3
    enableGradDOF = 1;
else
    enableGradDOF = 0;
end

%[node, elem, saved, forces, MPC, elementSelections] = DCB_2D(geometry,load_d); % Mesh the DCB specimen without cohesive elements 
[node, elem, saved, forces, MPC, elementSelections] = DCB_2D_CZ(geometry, load_d, enableGradDOF); % Mesh the DCB specimen using cohesive elements
[Dsave, saved_time, saved_forces, saved_displacements, history, DOF] = solver(node,elem,saved,forces,MPC,elementSelections,elementOptions,materialOptions,solverOptions);

%% Plotting
% Plot deformed result
plot_deformed_2d(node,elem,Dsave(:,end),DOF)

% Plot the cohesive elements
fig2 = figure;
elem_CZ = elem(elementSelections(:,1) == 2,:); % Get cohesive elements
D = Dsave(:,end); % The results for the last load step
plot_elem_opening(fig2,elem_CZ,node,D,DOF,elementOptions(2).elemID,'r');

% Plot force-displacement curve
fig3 = figure;
plot(saved_time,saved_forces(1,:),'s-');

% Plot interface tractions
fig4 = figure;
interface_elem = elem(elementSelections(:,1) == 2,:); % Get cohesive elements
interface_history = history(elementSelections(:,1) == 2,end); % Get history for the last load step
plot_interface_tractions(fig4, interface_elem, node, D, interface_history, DOF, material, elementOptions(2).elemID, 'r');