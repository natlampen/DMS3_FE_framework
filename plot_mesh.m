function plot_mesh(node, elem)
% This function plots the mesh from the inputs( nodes, elem)
    print_elem_no = 0;
    print_node_no = 1;
    
    n_nodes = length(node(:,1));
    n_elem = length(elem(:,1));
    fig1 = figure;
    hold on
    scatter(node(:,1),node(:,2),'r.')
    % Plot elements
    x_elem_text = zeros(n_elem,1);
    y_elem_text = zeros(n_elem,1);
    for j = 1:n_elem
        x_coords = node(elem(j,:),1);
        y_coords = node(elem(j,:),2);
        fill(x_coords,y_coords,'y');
        x_elem_text(j) = mean(x_coords);
        y_elem_text(j) = mean(y_coords);
    end
    % Print node numbers
    if print_node_no == 1
        for j = 1:n_nodes
            text(node(j,1),node(j,2),num2str(j),'fontsize',8);
        end
    end
    % Print element numbers
    if print_elem_no == 1
        for j = 1:n_elem
            text(x_elem_text(j),y_elem_text(j),num2str(j),'fontsize',4,'HorizontalAlignment','center');
        end
    end
    xlim([min(node(:,1))-20e-3 max(node(:,1))+20e-3]);
    ylim([min(node(:,2))-20e-3 max(node(:,2))+20e-3]);
    axis equal
end

