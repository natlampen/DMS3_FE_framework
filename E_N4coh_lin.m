function [k, R_in, history] = E_N4coh_lin(node, u, elemData, matModel, matData, history)
% This function formulates a cohesive zone interface element

K = matData(1);
Gc = matData(2);
tau_0 = matData(3);

% Calculate onset and final openings
Delta_0 = tau_0/K;
Delta_f = 2*Gc/tau_0;

% Remove all z-coords from input and reshape to column vector
node = reshape(node(:,1:2)',8,1);
u = reshape(u(:,1:2)',8,1);

% Allocate space for the element stiffness matrix and internal force vector
k = zeros(8);
R_in = zeros(8,1);

b = elemData(1);

% Setup integration points
if elemData(2) == 1
    sample_points = [0];
    gauss_weight = 2;
elseif elemData(2) == 2
    sample_points = [-0.577350269189626, 0.577350269189626];
    gauss_weight = [1, 1];
elseif elemData(2) == 3
    sample_points = [-0.774596669241483, 0.774596669241483, 1];
    gauss_weight = [0.555555555555556, 0.555555555555556, 0.888888888888889];
elseif elemData(2) == 20
    sample_points = [-0.9931285992,-0.9639719273,-0.9122344283,-0.8391169718,-0.7463319065, ...
                     -0.6360536808,-0.5108670020,-0.3737060887,-0.2277858511,-0.07652652113, ...
                      0.07652652113,0.2277858511,0.3737060887,0.5108670020,0.6360536808, ...
                      0.7463319065,0.8391169718,0.9122344283,0.9639719273,0.9931285992];
    gauss_weight = [0.01761400797,0.04060143030,0.06267204900,0.08327674140,0.1019301200, ...
                    0.1181945322,0.1316886387,0.1420961093,0.1491729865,0.1527533872, ...
                    0.1527533872,0.1491729865,0.1420961093,0.1316886387,0.1181945322, ...
                    0.1019301200,0.08327674140,0.06267204900,0.04060143030,0.01761400797];
end

% Setup the rotation matrix - This is constant due to linear shape
% functions
dN1 = -1/2;
dN2 = 1/2;
dN = [
    dN1 0   dN2 0  ;
    0   dN1 0   dN2];

g = 1/2*dN*(node(1:4) + node(5:8) + u(1:4) + u(5:8));
J = norm(g);
RotMat = 1/J*[g(1) -g(2); g(2) g(1)];

% Allocate memory for damage
damage = zeros(elemData(2),1);

if size(history) == size([])
    history = zeros(elemData(2),1);
end

% Use gauss quadrature for integration
n_sample_points = size(sample_points,2);
for ii = 1:n_sample_points
    xi = sample_points(ii);
    
    % Calculate local separations
    N1 = 1/2*(1-xi);
    N2 = 1/2*(1+xi);
    
    N = [
        N1 0  N2 0;
        0  N1 0 N2];
    N = [-N N];
    
    Delta = RotMat'*N*u;
    
    % Calculate damage for the integration point
    Delta_n = sqrt( (1/2*(abs(Delta(2))+Delta(2)))^2 + Delta(1)^2);
    Delta_t = Delta_0*Delta_f/(Delta_f - history(ii)*(Delta_f - Delta_0));
    if Delta_n < Delta_t
        damage(ii) = history(ii);
    elseif Delta_n < Delta_f
        damage(ii) = Delta_f/Delta_n*(Delta_n - Delta_0)/(Delta_f - Delta_0);
    else
        damage(ii) = 1;
    end
    
    % Calculate traction forces
    if Delta(2) >= 0
        tau = (1-damage(ii))*K*Delta;
    else
        tau = [(1-damage(ii))*K*Delta(1); K*Delta(2)];
    end
    
    Dtan = zeros(2);
    if abs(Delta(1)) < Delta_t
        Dtan(1,1) = (1-history(ii))*K;
    elseif abs(Delta(1)) < Delta_f
        Dtan(1,1) = -K*Delta_0/(Delta_f-Delta_0);
    else
        Dtan(1,1) = 0;
    end
    
    if Delta(2) < 0
        Dtan(2,2) = K;
    elseif Delta(2) < Delta_t
        Dtan(2,2) = (1-history(ii))*K;
    elseif Delta(2) < Delta_f
        Dtan(2,2) = -K*Delta_0/(Delta_f-Delta_0);
    else
        Dtan(2,2) = 0;
    end
     
    k = k + N'*RotMat*Dtan*RotMat'*N*J*gauss_weight(ii)*b;
    R_in = R_in + N'*RotMat*tau*J*gauss_weight(ii)*b;
end

% Save damage
history = damage;

end

