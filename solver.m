function [Dsave, saved_time, saved_forces, saved_displacements, history, DOF] = solver(node,elem,saved,forces,MPC,elementSelections,elementOptions,materialOptions,solver)
% Element formulations
availableElem = struct();

availableElem(1).handle = @E_N4iso_lin;
availableElem(1).nDOF = 2;
availableElem(1).DOF = {'UX','UY'};

availableElem(2).handle = @E_N4coh_lin;
availableElem(2).nDOF = 2;
availableElem(2).DOF = {'UX','UY'};

availableElem(3).handle = @E_N4coh_mpc;
availableElem(3).nDOF = 4;
availableElem(3).DOF = {'UX','UY','GRADX','GRADY'};

% Material models
availableMat = struct();

availableMat(1).handle = @M_LinearElastic;

% Determine node and element count
n_node = size(node,1);
n_elem = size(elem,1);

% Setup bookkeeping for tracking of indexes of the different DOF for nodes
clear DOF % Clear the DOF struct to prevent unintended behavior which will happen if DOF is defined
uniqueElements = unique([elementOptions(elementSelections(:,1)).elemID]); % Determine which elements are used in the model
uniqueDOFs = unique([availableElem(uniqueElements).DOF]); % Determine which types of DOFs are present in the model
DOF(n_node) = cell2struct(cell(size(uniqueDOFs)), uniqueDOFs, 2); % Create a struct array of size n_node x 1.
nextDOFindex = 1; % This is the first index to use
% Loop over all elements and for each node of that element add the DOFindex
% to the DOF struct if not already defined.
for ii = 1:n_elem
    elementNodes = elem(ii,:); % Determine the nodes belonging to element ii
    elementDOFs = availableElem(elementOptions(elementSelections(ii)).elemID).DOF; % Determine degress of freedom for element ii
    for jj = elementNodes
        for kk = elementDOFs
            if not(size(DOF(jj).(kk{1}))) % If the size of DOF kk for node jj is zero, the index is defined and is assigned nextDOFindex
                DOF(jj).(kk{1}) = nextDOFindex;
                nextDOFindex = nextDOFindex+1;
            end
        end
    end
end
n_DOF = nextDOFindex - 1; % Total amount of DOF of the model

% Setup the external force vector. This is constant with respect to
% displacements
R_ext = zeros(n_DOF,1); % Allocate space for the external force vector
for ii = forces % Loop over the applied forces and add to the external force vector
    R_ext(DOF(ii.node).(ii.DOF)) = ii.force;
end

% Allocate displacement vectors
D = zeros(n_DOF,1); % Displacement vector
Dsave = zeros(n_DOF,solver.max_iter+1); % Matrix for saving displacement vectors

% Setup the MPC matrix
n_MPC = size(MPC,2);
lambda = zeros(n_MPC,1);
lambda0 = zeros(n_MPC,1);
Q_ext = zeros(n_MPC,1);
C = zeros(n_MPC,n_DOF);
for ii = 1:n_MPC
    for jj = 1:size(MPC(ii).coef,2)
        C(ii,DOF(MPC(ii).node(jj)).(MPC(ii).DOF{jj})) = MPC(ii).coef(jj);
    end
    Q_ext(ii) = MPC(ii).constant;
end

% Setup data saving feature
saved_index = false(n_DOF,1); % index of the saved data
for ii = saved
    saved_index(DOF(ii.node).(ii.DOF)) = true;
end
n_saved = sum(saved_index);
saved_time = zeros(1,solver.max_iter);
saved_forces = zeros(n_saved,solver.max_iter);
saved_displacements = zeros(n_saved,solver.max_iter);

% Allow for element history. This cell is used to save the history of the
% last converged step.
history = cell(n_elem,solver.max_iter+1); % max_iter+1 to also have timestep = 0

% Start solver iterations
stepsize = solver.stepsize;
time = stepsize;
substep_iter = 1;
substep = 1;
for iter = 1:solver.max_iter
    % Reset history from last time step
    step_history = cell(n_elem,1);
    
    % Reset internal force vector
    R_in = zeros(n_DOF,1);
    
    % Allocate space for the tangent stiffness
    K = zeros(n_DOF);
    
    % Loop over all elements and add their local stiffness to the global
    % stiffness matrix
    n_elem = size(elem,1);
    for ii = 1:n_elem
        n_elnodes = size(elem(ii,:),2);
        elementDOFs = availableElem(elementOptions(elementSelections(ii)).elemID).DOF;
        n_elemDOFs = size(elementDOFs,2);
        el_index = zeros(n_elnodes,n_elemDOFs);
        for jj = 1:n_elemDOFs
            % el_index contains one row per node
            el_index(:,jj) = [DOF(elem(ii,:)).(elementDOFs{jj})]';
        end
        
        el_node = node(elem(ii,:),:); % Element nodal coordinates
        el_u = D(el_index);
        
        elementOptionsID = elementSelections(ii,1);
        materialOptionsID = elementSelections(ii,2);
        elemHandle = availableElem(elementOptions(elementOptionsID).elemID).handle; % The elemHandle represents the element formulation
        matModel = availableMat(materialOptions(materialOptionsID).matID).handle; % The matModel represents the material model
        
        % Evaluate the tangent stiffness of the element
        [k, el_R_in ,el_history] = feval(elemHandle, el_node, el_u, elementOptions(elementOptionsID).options, matModel, materialOptions(materialOptionsID).options, history{ii,substep});
        step_history{ii} = el_history;
        
        % Add element stiffness matrix to global stiffness matrix
        el_index = reshape(el_index',n_elnodes*n_elemDOFs,1);
        K(el_index,el_index) = K(el_index,el_index) + k;
        
        % Add element internal force vector to the global
        R_in(el_index) = R_in(el_index) + el_R_in;
    end
    % -------------------------------------------------------------------------
    
    % Calculate the residual
    R = R_in + C'*lambda - R_ext*time;
    Q = C*D - solver.perb*lambda - Q_ext*time;
    
    % Check tolerance of the current solution
    if norm([R; Q]) < solver.tol % The substep has converged
        fprintf('Time %f has converged in %d iteration(s)\n', time, substep_iter);

        % Automatic load stepping - Increase step if <= 4
        if substep_iter <= 4
            stepsize = min([stepsize*2, solver.max_stepsize]);
        end
        
        % Save load step data
        history(:,substep+1) = step_history; % Substep has converged - Update history
        Dsave(:,substep+1) = D;
        lambda0 = lambda;
        if exist('saved','var')
            saved_time(substep) = time;
            saved_forces(:,substep) = R_in(saved_index);
            saved_displacements(:,substep) = D(saved_index);
        end
        
        
        if time == 1
            fprintf('Solution done\n');
            history = history(:,2:substep+1); % History is offset by 1 to include timestep 0;
            Dsave = Dsave(:,2:substep+1); % Dsave is offset by 1 to include timestep 0;
            saved_time = saved_time(1:substep);
            saved_forces = saved_forces(1:substep);
            saved_displacements = saved_displacements(1:substep);
            break
        elseif time+stepsize >= 1
            stepsize = 1-time;
            time = 1;
        else
            time = time + stepsize;
        end
        
        substep_iter = 0;
        substep = substep + 1;
    end
    
    if substep_iter >= 12
        if stepsize == solver.min_stepsize
            fprintf('Could not converge using minimum stepsize\n');
            history = history(:,2:substep);
            Dsave = Dsave(:,2:substep);
            saved_time = saved_time(1:substep-1);
            saved_forces = saved_forces(1:substep-1);
            saved_displacements = saved_displacements(1:substep-1);
            break
        end
        time = time - stepsize;
        stepsize = max([stepsize/2, solver.min_stepsize]);
        time = time + stepsize;
        substep_iter = 0;
        D = Dsave(:,substep);
        lambda = lambda0;
    else
        % Update D
        delta_D = [K C'; C -solver.perb*eye(n_MPC)] \ -[R; Q];
        D = D + delta_D(1:n_DOF);
        lambda = lambda + delta_D(n_DOF+1:end);
    end
    substep_iter = substep_iter + 1;
end