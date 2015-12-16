function plot_elem_opening(fig, elem, node, D, DOF, elemID,color)
% Plot gausspoints
plot_gauss = 0;

% Number of evaluation points when plotting the deformation field
n_eval = 30;

% The integration points for 2nd order gauss integration. The integration
% order is not likely to change in this project
intPoints = [-0.577350269189626, 0.577350269189626];

% Remove z-coordinate
node = node(:,1:2);

% Make plot window active
figure(fig)
hold on

n_elem = size(elem,1);
for ii = 1:n_elem
    % Plot element nodes
    elem_nodes_x = node(elem(ii,:),1); % Select nodal coordinates in x related to element ii
    elem_nodes_y = node(elem(ii,:),2); % Select nodal coordinates in y related to element ii
    elem_nodes_ux = D([DOF(elem(ii,:)).UX]); % Select nodal displacements in x related to element ii
    elem_nodes_uy = D([DOF(elem(ii,:)).UY]); % Select nodal displacements in y related to element ii
    if elemID == 3
        elem_nodes_gradx = D([DOF(elem(ii,:)).GRADX]); % Select nodal displacements in x related to element ii
        elem_nodes_grady = D([DOF(elem(ii,:)).GRADY]); % Select nodal displacements in y related to element ii
    end
    scatter(elem_nodes_x+elem_nodes_ux,elem_nodes_y+elem_nodes_uy,[],color); % Plot nodes
        
    %
    % Plot integrations points
    %
    L = elem_nodes_x(2) - elem_nodes_x(1); % Length of midplane - Assuming no/minimal deformation
    n_intPoints = length(intPoints);
    if plot_gauss
        for jj = 1:n_intPoints
            xi = intPoints(jj);
            
            Ng1 = (1-xi)/2;
            Ng2 = (1+xi)/2;
            Ng = [Ng1 0 Ng2 0; 0 Ng1 0 Ng2]; % Evaluate the geometric shapefunction
            c = reshape([elem_nodes_x elem_nodes_y]',8,1); % Format the nodal coordinates appropriate for the calculation of integration point coordinates
            coord_geom = [Ng Ng]*c/2;
            
            if elemID == 3
                % Evaluate the shapefunctions for the MPC approach
                x = (1+xi)*L/2;
                Nd1 = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
                Nd2 = x - (2*x^2)/L + x^3/L^2;
                Nd3 = (3*x^2)/L^2 - (2*x^3)/L^3;
                Nd4 = x^3/L^2 - x^2/L;
                
                % Evaluate the nodal displacements
                d = reshape([elem_nodes_ux elem_nodes_uy elem_nodes_gradx elem_nodes_grady]',16,1);
                Nd = [Nd1 0 Nd2 0 Nd3 0 Nd4 0; 0 Nd1 0 Nd2 0 Nd3 0 Nd4];
                coord_disp = [Nd Nd]*d/2;
            elseif elemID == 3
                % Evaluate nodal displacements
                d = reshape([elem_nodes_ux elem_nodes_uy]',8,1); % Format the displacements appropriate for the calculate of the integration point coordinates
                coord_disp = [Ng Ng]*d/2;
            else
                disp('Undefined elemID');
                return
            end
            
            % Calculate gauss point coordinate
            intPointCoord = coord_geom + coord_disp;
            plot([intPointCoord(1) intPointCoord(1)],[-1 1],'--k');
        end
    end
    
    %
    % Plot displacement field between nodes
    %
    coords_lower = zeros(2,n_eval);
    coords_upper = zeros(2,n_eval);
    xiList = linspace(-1,1,n_eval);
    for jj = 1:n_eval
        xi = xiList(jj);
        
        Ng1 = (1-xi)/2;
        Ng2 = (1+xi)/2;
        Ng = [Ng1 0 Ng2 0; 0 Ng1 0 Ng2]; % Evaluate the geometric shapefunction
        c = reshape([elem_nodes_x elem_nodes_y]',8,1);
        
        coords_lower(:,jj) = Ng*c(1:4);
        coords_upper(:,jj) = Ng*c(5:8);
        
        if elemID == 3
            x = (1+xi)*L/2;
            Nd1 = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
            Nd2 = x - (2*x^2)/L + x^3/L^2;
            Nd3 = (3*x^2)/L^2 - (2*x^3)/L^3;
            Nd4 = x^3/L^2 - x^2/L;
            Nd = [Nd1 0 Nd2 0 Nd3 0 Nd4 0; 0 Nd1 0 Nd2 0 Nd3 0 Nd4];
            d = reshape([elem_nodes_ux elem_nodes_uy elem_nodes_gradx elem_nodes_grady]',16,1);
            
            coords_lower(:,jj) = coords_lower(:,jj) + Nd*d(1:8);
            coords_upper(:,jj) = coords_upper(:,jj) + Nd*d(9:16);
        elseif elemID == 2
            d = reshape([elem_nodes_ux elem_nodes_uy]',8,1);
            
            coords_lower(:,jj) = coords_lower(:,jj) + Ng*d(1:4);
            coords_upper(:,jj) = coords_upper(:,jj) + Ng*d(5:8);
        else
            disp('Undefined elemID');
            return
        end
    end
    
    plot(coords_lower(1,:),coords_lower(2,:),'color',color);
    plot(coords_upper(1,:),coords_upper(2,:),'color',color);
end