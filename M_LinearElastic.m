function C = M_LinearElastic(matData)
% This function returns the constitutive matrix for a linear elastic
% material.
% For isotropic
% MatData = [1, Youngs Modulus, Poisson Ratio]

if matData(1) == 1 % Isotropic
    E = matData(2);
    nu = matData(3);
    
    C11 = E/(1+nu)*((1-nu)/(1-2*nu));
    C12 = E/(1+nu)*(nu/(1-2*nu));
    C44 = (C11-C12)/2;
    
    C = [
        C11 C12 C12 0   0   0;
        C12 C11 C12 0   0   0;
        C12 C12 C11 0   0   0;
        0   0   0   C44 0   0;
        0   0   0   0   C44 0;
        0   0   0   0   0   C44];
else
    disp('Material Model Error');
end
    
end

