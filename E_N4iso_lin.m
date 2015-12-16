function [k, R_in, history] = E_N4iso_lin(node, u, elemData, matModel, matData, el_history)
% This function returns the local stiffness matrix of a Q4-isoparametric
% structural element based on the nodal coordinates and element parameters.
% Element options: (PLANE-STRESS)
%  1 - Element thickness
%  2 - Integration order
%  3 - Stress assumption (1 = Plane-Stress)

% History is not used for this element
history = [];

% Allocate space for the element stiffness matrix
k = zeros(8);

% Reshape displacements vector
u = reshape(u(:,1:2)',8,1);

% Setup integration points
if elemData(2) == 1 % Reduced integration
    sample_points = [0];
    gauss_weight = 2;
elseif elemData(2) == 2 % Full integration
    sample_points = [-1/sqrt(3), 1/sqrt(3)];
    gauss_weight = 1;
end

% Setup the constitutive matrix
E = feval(matModel, matData);
if elemData(3) == 1 % Planestress
    E11 = E(1,1) - E(1,2)^2/E(1,1);
    E12 = E(1,2) - E(1,2)^2/E(1,1);
    E33 = E(4,4);
    E = [
        E11 E12 0;
        E12 E11 0;
        0   0   E33];
end
        
n_sample_points = size(sample_points,2);
for ii = 1:n_sample_points
    eta = sample_points(ii);
    for jj = 1:n_sample_points
        xi = sample_points(jj);
        
        shape_diff = 1/4*[-(1-eta), (1-eta), (1+eta), -(1+eta); -(1-xi), -(1+xi), (1+xi), (1-xi)];
        exp_shape_diff = zeros(4,8);
        exp_shape_diff(1:2,1:2:8) = shape_diff;
        exp_shape_diff(3:4,2:2:8) = shape_diff;
        
        J = shape_diff*node(:,1:2);
        detJ = (J(1,1)*J(2,2)-J(2,1)*J(1,2));
        Gamma = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
        exp_Gamma = zeros(4);
        exp_Gamma(1:2,1:2) = Gamma;
        exp_Gamma(3:4,3:4) = Gamma;
        
        B = [1 0 0 0;0 0 0 1; 0 1 1 0]*exp_Gamma*exp_shape_diff;
        k = k + (B')*E*B*elemData(1)*detJ*gauss_weight^2;
    end
end

% Evaluate the internal force vector for the element
R_in = k*u;

end

