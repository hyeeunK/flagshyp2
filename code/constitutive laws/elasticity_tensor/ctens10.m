%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (in Voigt notation) for material type 1.
%--------------------------------------------------------------------------
function c = ctens10(kinematics,properties,cons)
mu1         = properties(2);
mu2         = properties(2);
lambda     = properties(3);
k = lambda + 2*mu1/3;
J          = kinematics.J;

C = kinematics.F' * kinematics.F;
C_bar = J^(-2/3) * C;
I1 = trace(C);
I1_bar = J^(-2/3) * I1;
I2 = 1/2 * (I1^2 - trace(C^2));
I2_bar = J^(-4/3) * I2;

p = cons.I - 1/3 * C_bar^(-1)*C_bar^(-1)';
p_tilde = k * (2*J -1);
S_bar = (mu1 + mu2*I1_bar)*cons.I - mu2*C_bar;

% compute pk2
S_vol = J^(1/3) * (J-1) * k * C_bar^(-1);
S_iso = J^(-2/3) * (-1/3*(mu1*I1_bar + 2*mu2*I2_bar)*C_bar^(-1) + (mu1 + mu2*I1_bar)*cons.I - mu2*C_bar);
S = S_vol + S_iso;

c_vol = J^(-1/3)*(2*J-1)*k*C_bar^(-1)*C_bar^(-1)' - 2*J^(-1/3)*(J-1)*k*C_bar;
c_iso = p*C_bar*p' + 2/3*trace(J^(-2/3)*S_bar)*p_tilde -2/3*J^(-2/3)*(C_bar^(-1)*S_iso' + S_iso*C_bar^(-1)');
c = c_vol + c_iso;

end


