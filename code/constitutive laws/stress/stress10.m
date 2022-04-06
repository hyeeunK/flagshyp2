%--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor for material type 1.
%--------------------------------------------------------------------------
function Cauchy = stress10(kinematics,properties,cons)
mu1              = properties(2);
mu2              = properties(2);
lambda          = properties(3);
J               = kinematics.J;
b               = kinematics.b;
k = lambda + 2*mu1/3;
E = 1/2*(kinematics.F' * kinematics.F - cons.I);

C = kinematics.F' * kinematics.F;
C_bar = J^(-2/3) * C;
I1 = trace(C);
I1_bar = J^(-2/3) * I1;
I2 = 1/2 * (I1^2 - trace(C^2));
I2_bar = J^(-4/3) * I2;

% compute pk2
S_vol = J^(1/3) * (J-1) * k * C_bar^(-1);
S_iso = J^(-2/3) * (-1/3*(mu1*I1_bar + 2*mu2*I2_bar)*C_bar^(-1) + (mu1 + mu2*I1_bar)*cons.I - mu2*C_bar);
S = S_vol + S_iso;

% push forward
Cauchy          = 1/J * kinematics.F * S * kinematics.F';
end