function plot_interface_tractions(fig, interface_elem, node, D, interface_history, DOF, material, elemID, color)
figure(fig); hold on; % Selection figure for plotting

% Get number of interface elements
n_elem = size(interface_elem,1);

% Number of gauss points
n_gp = size(interface_history{1},1);

% Define gauss point coordinates
if n_gp == 2
    gp_coord = [-0.577350269189626, 0.577350269189626];
elseif n_gp == 20
    gp_coord = [-0.9931285992,-0.9639719273,-0.9122344283,-0.8391169718,-0.7463319065, ...
        -0.6360536808,-0.5108670020,-0.3737060887,-0.2277858511,-0.07652652113, ...
        0.07652652113,0.2277858511,0.3737060887,0.5108670020,0.6360536808, ...
        0.7463319065,0.8391169718,0.9122344283,0.9639719273,0.9931285992];
end

tractions = zeros(n_gp*n_elem, 1); % Allocation space for the tractions - It is calculated in the gauss points
x_coord = zeros(n_gp*n_elem, 1); % Allocate space for the x-coordinate associated with the tractions
for ii = 1:n_elem
    element_nodes = interface_elem(ii,:); % Get nodes for interface element ii
    element_history = interface_history{ii};
    
    % Get nodal coordinates
    nodal_xx = node(element_nodes,1);
    
    % Get displacements of element nodes
    nodal_ux = D([DOF(element_nodes).UX]);
    nodal_uy = D([DOF(element_nodes).UY]);
    if elemID == 3
        nodal_gradx = D([DOF(element_nodes).GRADX]);
        nodal_grady = D([DOF(element_nodes).GRADY]);
        u = reshape([nodal_ux nodal_uy nodal_gradx nodal_grady]',16,1);
    else
        u = reshape([nodal_ux nodal_uy]',8,1);
    end
    
    for jj = 1:n_gp
        % Calculate the shapefunctions
        xi = gp_coord(jj);
        
        if elemID == 3
            L = nodal_xx(4) - nodal_xx(3);
            x = (1+xi)*L/2;
            
            % Beam shapefunctions
            N1 = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
            N2 = x - (2*x^2)/L + x^3/L^2;
            N3 = (3*x^2)/L^2 - (2*x^3)/L^3;
            N4 = x^3/L^2 - x^2/L;
            N = [N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4];
        else
            N1 = (1-xi)/2;
            N2 = (1+xi)/2;
            N = [N1 0 N2 0; 0 N1 0 N2];
        end
        
        % Because of DCB the opening displacements can be calculated as:
        % (The rotation matrix becomes the identity matrix)
        Delta = [-N N]*u;
        
        % Calculate tractions
        gp_traction = (1-element_history(jj))*material.penalty_stiffness*Delta - element_history(jj)*material.penalty_stiffness*(abs(Delta(2)) - Delta(2))/2;
        tractions((ii-1)*n_gp+jj) = gp_traction(2); % Only the mode I component
        
        % Geometry shapefunctions
        N1 = (1-xi)/2;
        N2 = (1+xi)/2;
        
        % Calculate the x-coordinate
        x_coord((ii-1)*n_gp+jj) = [N1 N2 N1 N2]*nodal_xx/2;
    end
end

plot(x_coord,tractions,'color',color); %,'color',plotColors(1,:));
end