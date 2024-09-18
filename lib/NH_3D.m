function [T,R,Lt,T_vol] = NH_3D(F,params,eq,J,X,params_F)

mu = X*params(1);
kappa = X*params(2);
nu = X*params(3);
temp   = params_F(1);
temp0  = params_F(2);
alphaT = params_F(3);
zita   = params_F(4);
alphaZ = params_F(5);

% Jt     = 1 + alphaT*(temp - temp0);
% Jz     = 1 + alphaZ*zita;
Jt = 1.0;
Jz = 1.0;
Jm     = J/(Jt*Jz);
%  Fstar = J^(-1/3) * F;
Fstar = F;
bstar = Fstar * Fstar';
devbstar = Dev(bstar);
Ibar1 = trace(bstar);
lc    = sqrt(Ibar1 / 3);
Lt     = sqrt(X*(lc^2 - 1) + 1);
I         = eye(3);
C = F'*F;

if eq == 1
    T = mu/J * devbstar + 0.5*kappa*(Jm - 1/Jm)/(Jt*Jz)*I;
    T_vol = 0.5*kappa*(Jm - 1/Jm)/(Jt*Jz)*I;
else
    T = mu/J * devbstar;
end

V           = zeros(3);
lnV         = zeros(3);
be          = F*F';


if eq == 1
    [Q,lambda2] = eig(be);
    lambda      = sqrt(diag((lambda2)));
    for i = 1:3
        ni  = Q(:,i); 
        V   = V + lambda(i) * (ni*ni');    
        lnV = lnV + log(lambda(i)) * (ni*ni');
    end
    R = (I/V)*F;
else
    R = eye(3);
end


end

