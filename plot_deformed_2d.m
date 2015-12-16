function plot_deformed_2d(node, elem, D, DOF)
% This function plots the mesh from the inputs( nodes, elem)
    print_node_no = 0;
    print_elem_no = 0;

    n_nodes = length(node(:,1));
    n_elem = length(elem(:,1));
    deformed_nodes = node(:,1:2) + D([DOF.UX;DOF.UY]');
    fig1 = figure;
    hold on
    % Plot elements
    x_elem_text = zeros(n_elem,1);
    y_elem_text = zeros(n_elem,1);
    for j = 1:n_elem
        x_coords = deformed_nodes(elem(j,:),1);
        y_coords = deformed_nodes(elem(j,:),2);
        fill(x_coords,y_coords,'y');
        x_elem_text(j) = mean(x_coords);
        y_elem_text(j) = mean(y_coords);
    end
    % Plot nodes
    scatter(deformed_nodes(:,1),deformed_nodes(:,2),'r.');
    % Print node numbers
    if print_node_no == 1
        for j = 1:n_nodes
            text(deformed_nodes(j,1),deformed_nodes(j,2),num2str(j),'fontsize',8);
        end
    end
    % Print element numbers
    if print_elem_no == 1
        for j = 1:n_elem
            text(x_elem_text(j),y_elem_text(j),num2str(j),'fontsize',4,'HorizontalAlignment','center');
        end
    end
    xlim([min(deformed_nodes(:,1))-20e-3 max(deformed_nodes(:,1))+20e-3]);
    ylim([min(deformed_nodes(:,2))-20e-3 max(deformed_nodes(:,2))+20e-3]);
    axis equal
end

