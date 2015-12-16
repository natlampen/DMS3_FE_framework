function [node, elem, saved, forces, MPC, elementSelections] = DCB_2D(geometry, load_d)

% Define mesh
DCB_length = geometry.DCB_length;
DCB_height = geometry.DCB_height;
DCB_width = geometry.DCB_width;
crack_length = geometry.DCB_a0;
DCB_elem_size_height = geometry.DCB_elem_size_height;
DCB_elem_size_length = geometry.DCB_elem_size_length;

%The first node is placed at (x,y) = (0,0)
x0 = 0;
y0 = 0;
z0 = 0;

% The body of the DCB specimen is created first.
body_length = DCB_length-crack_length;
n_elem_body_length = round(body_length / DCB_elem_size_length);
n_elem_body_height = round(DCB_height / DCB_elem_size_height);

% n_elem + 1 nodes are needed
node = [];
for ii = 1:n_elem_body_length+1
    for jj = 1:n_elem_body_height+1
        % This is not a very efficient way of during it since MATLAB has to reallocate space for every iteration
        node(end+1,:) = [x0 + (ii-1)*DCB_elem_size_length, y0 + (jj-1)*DCB_elem_size_height, z0];    
    end
end
elem = []; % This is the 4 node structural element
for ii = 1:n_elem_body_length
    for jj = 1:n_elem_body_height
        % Not very efficient
        node1 = ii + (jj-1) + (ii-1)*n_elem_body_height;
        elem(end+1,:) = [node1, node1+n_elem_body_height+1, node1+n_elem_body_height+2, node1+1];
    end
end

% The lower arm is created next
n_offs = length(node(:,1));
n_elem_arm_length = crack_length / DCB_elem_size_length;
n_elem_arm_height = DCB_height / 2 / DCB_elem_size_height;
for ii = 1:n_elem_arm_length
    for jj = 1:n_elem_arm_height+1
        % Not very effective
        node(end+1,:) = [x0 + body_length + (ii)*DCB_elem_size_length, y0 + (jj-1)*DCB_elem_size_height, z0];
    end
end
for ii = 1:n_elem_arm_length-1
    for jj = 1:n_elem_arm_height
        % Not very effective
        node1 = n_offs + ii + (jj-1) + (ii-1)*n_elem_arm_height;
        elem(end+1,:) = [node1, node1+n_elem_arm_height+1, node1+n_elem_arm_height+2, node1+1];
    end
end
for ii = 1:n_elem_arm_height
    node1 = (n_elem_body_length)*(n_elem_body_height+1)+ii;
    node2 = (n_elem_body_length+1)*(n_elem_body_height+1)+ii;
    elem(end+1,:) = [node1, node2, node2+1, node1+1];
end

% The upper arm is created
n_offs = length(node(:,1));
for ii = 1:n_elem_arm_length
    for jj = 1:n_elem_arm_height+1
        % Not very effective
        node(end+1,:) = [x0 + body_length + (ii)*DCB_elem_size_length, y0 + DCB_height/2 + (jj-1)*DCB_elem_size_height, z0];
    end
end
for ii = 1:n_elem_arm_length-1
    for jj = 1:n_elem_arm_height
        % Not very effective
        node1 = n_offs + ii + (jj-1) + (ii-1)*n_elem_arm_height;
        elem(end+1,:) = [node1, node1+n_elem_arm_height+1, node1+n_elem_arm_height+2, node1+1];
    end
end
for ii = 1:n_elem_arm_height
    node1 = (n_elem_body_length)*(n_elem_body_height+1)+n_elem_body_height/2+ii;
    node2 = n_offs+ii;
    elem(end+1,:) = [node1, node2, node2+1, node1+1];
end

% Create structs for forces and MPCs
forces = struct('node',{},'DOF',{},'force',{});
MPC = struct('node',{},'DOF',{},'coef',{},'constant',{});

% Locate nodes where forces should be applied
upper_node = find(abs(node(:,1)-DCB_length)<1e-5 & abs(node(:,2)-DCB_height)<1e-5);
lower_node = find(abs(node(:,1)-DCB_length)<1e-5 & abs(node(:,2)-0<1e-5));

% Apply constraints and forces
MPC(1).node = lower_node;
MPC(1).DOF = {'UX'};
MPC(1).coef = 1;
MPC(1).constant = 0;

MPC(2).node = lower_node;
MPC(2).DOF = {'UY'};
MPC(2).coef = 1;
MPC(2).constant = -load_d/2;

MPC(3).node = upper_node;
MPC(3).DOF = {'UX'};
MPC(3).coef = 1;
MPC(3).constant = 0;

MPC(4).node = upper_node;
MPC(4).DOF = {'UY'};
MPC(4).coef = 1;
MPC(4).constant = load_d/2;

% What to save?
saved = struct('node',{},'DOF',{});
saved(1).node = upper_node;
saved(1).DOF = 'UY';

elementSelections = ones(size(elem,1),2); % Select Q4isoparametric and linear elastic for specimen
