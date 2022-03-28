%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (in Voigt notation) for material type 1.
%--------------------------------------------------------------------------
function c = ctens9(kinematics,properties,cons)
mu         = properties(2);
lambda     = properties(3);
J          = kinematics.J;
c          = lambda*cons.IDENTITY_TENSORS.c1 + mu*cons.IDENTITY_TENSORS.c2;
end


