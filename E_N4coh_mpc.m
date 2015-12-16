function [k, R_in, history] = E_N4coh_mpc(node, u, elemData, matModel, matData, history)
% This function formulates a cohesive zone interface element

K = matData(1);
Gc = matData(2);
tau_0 = matData(3);

% Calculate onset and final openings
Delta_0 = tau_0/K;
Delta_f = 2*Gc/(tau_0);

% Allocate space for the element stiffness matrix and internal force vector
nDOF = 16;
k = zeros(nDOF);
R_in = zeros(nDOF,1);

b = elemData(1);

% Reshape u to a column vector
u = reshape(u',nDOF,1);

% Reshape node to a column vector disregarding z-coord
node = reshape(node(:,1:2)',8,1);

% Calculate element length - Should this be the length of the deformed
% midplane - Like the rotation matrix is for the deformed geometry?
L = sqrt((node(7)-node(5))^2 + (node(8)-node(6))^2);

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

% Allocate memory for damage
damage = zeros(elemData(2),1);

if size(history) == size([])
    history = zeros(elemData(2),1);
end

% Use gauss quadrature for integration
n_sample_points = size(sample_points,2);
for ii = 1:n_sample_points
    xi = sample_points(ii);
    x = 1/2*(1+xi)*L;
    
    % Beam shapefunctions and derivatives with respect to x
    N1 = (2*x^3)/L^3 - (3*x^2)/L^2 + 1;
    dN1 = 6*x^2/L^3 - 6*x/L^2 + 1;
    N2 = x - (2*x^2)/L + x^3/L^2;
    dN2 = 1 - 4*x/L + 3*x^2/L^2;
    N3 = (3*x^2)/L^2 - (2*x^3)/L^3;
    dN3 = 6*x/L^2 - 6*x^2/L^3;
    N4 = x^3/L^2 - x^2/L;
    dN4 = 3*x^2/L^2 - 2*x/L;
    
    % Linear shapefunctions and derivatives with respect to x
    %N1lin = 1 - x/L;
    dN1lin = -1/L;
    %N2lin = x/L;
    dN2lin = 1/L;
    
    % Calculate the jacobian
    dN = [dN1 0 dN2 0 dN3 0 dN4 0; 0 dN1 0 dN2 0 dN3 0 dN4];
    dNlin = [dN1lin 0 dN2lin 0; 0 dN1lin 0 dN2lin];
    %J = norm( 1/2*([dN dN]*L/2*u + [dNlin dNlin]*L/2*node) );
    J = L/2;

    % Setup the rotation matrix - This is okay for DCB specimen
    RotMat = eye(2);
    
    % Calculate local separations
    N_pm = [N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4];
    N = [-N_pm N_pm];
    
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
     
    k = k + N'*RotMat*Dtan*RotMat'*N*J*b*gauss_weight(ii);
    R_in = R_in + N'*RotMat*tau*J*b*gauss_weight(ii);
end

% Save damage
history = damage;

end

