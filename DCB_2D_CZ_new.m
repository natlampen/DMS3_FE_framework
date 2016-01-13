function [node, elem, saved, forces, MPC, elementSelections] = DCB_2D_CZ(geometry, load_d, enableGradDOF)
% Analysis parameters
% Geometric parameters
DCB_height = geometry.DCB_height;
DCB_length = geometry.DCB_length;
DCB_width = geometry.DCB_width;
DCB_elem_size_height = geometry.DCB_elem_size_height;
DCB_elem_size_length = geometry.DCB_elem_size_length;
DCB_a0 = geometry.DCB_a0; %DCB_elem_size_length*20;

%The first node is placed at (x,y) = (0,0)
x0 = 0;
y0 = 0;
z0 = 0;

% The lower part of the specimen is created first
n_elem_length = round(DCB_length / DCB_elem_size_length);
n_elem_height = round(DCB_height/ 2 / DCB_elem_size_height);
% n_elem + 1 nodes are needed
node = [];
for ii = 1:n_elem_length+1
    for jj = 1:n_elem_height+1
        % This is not a very efficient way of during it since MATLAB has to reallocate space for every iteration
        node(end+1,:) = [x0 + (ii-1)*DCB_elem_size_length, y0 + (jj-1)*DCB_elem_size_height, z0];    
    end
end
elem = []; % This is the 4 node structural element
for ii = 1:n_elem_length
    for jj = 1:n_elem_height
        % Not very efficient
        node1 = ii + (jj-1) + (ii-1)*n_elem_height;
        elem(end+1,:) = [node1, node1+n_elem_height+1, node1+n_elem_height+2, node1+1];
    end
end

% The upper part is created next
y0 = DCB_height/2;
n_offs = length(node(:,1));
for ii = 1:n_elem_length+1
    for jj = 1:n_elem_height+1
        % Not very effective
        node(end+1,:) = [x0 + (ii-1)*DCB_elem_size_length, y0 + (jj-1)*DCB_elem_size_height, z0];
    end
end
for ii = 1:n_elem_length
    for jj = 1:n_elem_height
        % Not very effective
        node1 = n_offs + ii + (jj-1) + (ii-1)*n_elem_height;
        elem(end+1,:) = [node1, node1+n_elem_height+1, node1+n_elem_height+2, node1+1];
    end
end

% Element selection
elementSelections = ones(size(elem,1),2); % Select Q4isoparametric and linear elastic for specimen

% Add cohesive elements between upper and lower part of the specimen
for ii = 1:(n_elem_length-DCB_a0/DCB_elem_size_length)
    % Not very effective assignment
    node1 = ii*(n_elem_height+1);
    node2 = (n_elem_height+1)*(n_elem_length+1)+ 1 + (ii-1)*(n_elem_height+1);
    elem(end+1,:) = [node1, node1+n_elem_height+1, node2, node2+n_elem_height+1];
end
elementSelections = [elementSelections; 2*ones(n_elem_length-DCB_a0/DCB_elem_size_length,2)]; % Select cohesive zone


% Create structs to contain applied forces and MPCs
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

% Curvature in nodes
if enableGradDOF
    node1 = (n_elem_height+1);
    node2 = 2*(n_elem_height+1);
    node3 = 3*(n_elem_height+1);
   
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UY','UY','UY','GRADY'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UX','UX','UX','GRADX'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    
    node1 = (n_elem_height+1)*(n_elem_length+1) + 1;
    node2 = (n_elem_height+1)*(n_elem_length+1) + 1 + (n_elem_height+1);
    node3 = (n_elem_height+1)*(n_elem_length+1) + 1 + 2*(n_elem_height+1);
    
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UY','UY','UY','GRADY'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UX','UX','UX','GRADX'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    
    for ii = 2:(n_elem_length-DCB_a0/DCB_elem_size_length)
        node1 = (ii-1)*(n_elem_height+1);
        node2 = ii*(n_elem_height+1);
        node3 = (ii+1)*(n_elem_height+1);
        
        MPC(end+1) = cell2struct({[node1,node1,node2,node3,node3],{'GRADY','UY','GRADY','UY','GRADY'},[2*DCB_elem_size_length,6,8*DCB_elem_size_length,-6,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
        MPC(end+1) = cell2struct({[node1,node1,node2,node3,node3],{'GRADX','UX','GRADX','UX','GRADX'},[2*DCB_elem_size_length,6,8*DCB_elem_size_length,-6,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
        
        node1 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii-2)*(n_elem_height+1);
        node2 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii-1)*(n_elem_height+1);
        node3 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii)*(n_elem_height+1);
    
        MPC(end+1) = cell2struct({[node1,node1,node2,node3,node3],{'GRADY','UY','GRADY','UY','GRADY'},[2*DCB_elem_size_length,6,8*DCB_elem_size_length,-6,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
        MPC(end+1) = cell2struct({[node1,node1,node2,node3,node3],{'GRADX','UX','GRADX','UX','GRADX'},[2*DCB_elem_size_length,6,8*DCB_elem_size_length,-6,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    end
    
    ii = n_elem_length-DCB_a0/DCB_elem_size_length+1;
    
    node1 = (ii-2)*(n_elem_height+1);
    node2 = (ii-1)*(n_elem_height+1);
    node3 = ii*(n_elem_height+1);
    
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UY','UY','UY','GRADY'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UX','UX','UX','GRADX'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    
    node1 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii-3)*(n_elem_height+1);
    node2 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii-2)*(n_elem_height+1);
    node3 = (n_elem_height+1)*(n_elem_length+1) + 1 + (ii-1)*(n_elem_height+1);
    
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UY','UY','UY','GRADY'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
    MPC(end+1) = cell2struct({[node1,node2,node3,node3],{'UX','UX','UX','GRADX'},[-1,4,-3,2*DCB_elem_size_length],0},{'node','DOF','coef','constant'},2);
end

saved = struct('node',{},'DOF',{});
saved(1).node = upper_node;
saved(1).DOF = 'UY';