%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (in Voigt notation) for material type 1.
%--------------------------------------------------------------------------
function c = ctens10(kinematics,properties,cons)
mu1         = properties(2);
mu2         = properties(3);
k           = properties(4);
J           = kinematics.J;

C = kinematics.F' * kinematics.F;
C_bar = J^(-2/3) * C;
I1 = trace(C);
I1_bar = J^(-2/3) * I1;
I2 = 1/2 * (I1^2 - trace(C^2));
I2_bar = J^(-4/3) * I2;

dimension = size(C_bar,1);

B_bar = C_bar';
CONSTANT.I = eye(dimension); 

term1 =  zeros(dimension,dimension,dimension,dimension);
term2 =  zeros(dimension,dimension,dimension,dimension);
term3 =  zeros(dimension,dimension,dimension,dimension);

for i = 1:dimension
    for j = 1:dimension
        for k = 1:dimension
            for l = 1:dimension
                term1(i,j,k,l) = term1(i,j,k,l) + B_bar(i,j)*B_bar(k,l)...
                                 - 1/2*(B_bar(i,k)*B_bar(j,l) + B_bar(i,l)*B_bar(j,k));
                term2(i,j,k,l) = term2(i,j,k,l) + (B_bar(i,j)*CONSTANT.I(k,l) + B_bar(k,l)*CONSTANT.I(i,j));
                term3(i,j,k,l) = term3(i,j,k,l) + (B_bar(i,j)^2*CONSTANT.I(k,l) + B_bar(k,l)^2*CONSTANT.I(i,j));
            end
        end
    end
end

c_vol = k*J*(2*J-1)*cons.IDENTITY_TENSORS.c1 - k*J*(J-1)*cons.IDENTITY_TENSORS.c1; % Eq.(A.5c)
c_iso = 2*mu2*(term1)...
        -2/3*(mu1+2*mu2*I1_bar)*(term2)...
        +4/3*mu2*(term3) + 2/9*(mu1*I1_bar + 4*mu2*I2_bar)*cons.IDENTITY_TENSORS.c1...
        +1/3*(mu1*I1_bar + 2*mu2*I2_bar)*cons.IDENTITY_TENSORS.c2;      % Eq.(A.7)

c = c_vol + c_iso;
end


