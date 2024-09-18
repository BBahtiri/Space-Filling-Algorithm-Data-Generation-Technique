function [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tauHat,gamDot,strain_now] = Fi_backward(F,F_b_t,F_b_t_dt,F_p_t,F_p_t_dt...
    ,F_ve_t,TV,R,params,d,TA_t,dt,de,J,X,params_F,strain)
global eps0
%mu1 = params(1); mu2 = params(2);
%C = params(3); m = params(4); Ad = params(5); tau0 = params(6);
%d0s = params(7); m_tau = params(8); a = params(9); b = params(10);
%sigma0 = params(11);

mu1 = params(1); mu2 = params(2); nu1 = params(3); nu2 = 0;
m = params(4); gamma_dot_0 = params(5); dG = params(6);
Ad = params(7); tau0 = params(8); d0s = params(9); m_tau = params(10);
a = params(11); b = params(12); sigma0 = params(13);
y0 = params(14); x0 = params(15); a_t = params(16); b_t = params(17);
% kappa1 = 500;
% kappa2 = 500;
% kappa1 = 6 * mu1;
% kappa2 = 6 * mu2;
kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
kappa2 = 0;
% a = 0;
% b = 0;
% sigma0 = 100000;
% xi = 0.01;
% kappa1 = 500;
% kappa2 = 500;
% kappa1 = 7 * mu1;
% kappa2 = 7 * mu2;
params_A = [mu1, kappa1,nu1];
params_B = [mu2, kappa2,nu2];
k = 1.380649e-23;
T = params_F(1);


% Fve = F*(F_p_t)^(-1);
% Fe = Fve*(F_b_t)^(-1);



T_hat_v     = R'*TV*R;
tau         = norm(Dev(T_hat_v),'fro');
devStress = Dev(T_hat_v);
% devTV = Dev(TV);
% tau = norm(devTV,'fro');
% lamCh = sqrt((trace(F_b_t*F_b_t')./3));
% lamFac = lamCh - 1.0 + xi;
% tauHat = ((1-d) / (d0s))^(m_tau) * tau0;
% tauHat = tau0 - ((d / d0s)^m_tau);
tauHat = y0+a_t./(1+exp(-(d-x0)./b_t));
% tauHat = tau0;
gamDot = gamma_dot_0*exp(dG/(k*T)*((tau/tauHat)^(m) - 1));

%gamDot = lamFac^C * (tau / (tauHat * 1))^m;
prefac = 0.0;

if tau ~= 0
    prefac = gamDot / tau;
    F_b_t_new = (prefac *devStress*F_b_t)*dt + F_b_t_dt;
else
    F_b_t_new = F_b_t;
end

V           = zeros(3);
lnV         = zeros(3);
be          = F*F';
[Q,lambda2] = eig(be);
lambda      = sqrt(diag((lambda2)));
for i = 1:3
    ni  = Q(:,i); 
    V   = V + lambda(i) * (ni*ni');    
    lnV = lnV + log(lambda(i)) * (ni*ni');
end
strain_now = norm(lnV,'fro');
% strain1 = strain_now + strain;
strainDot = abs((strain_now - strain) / dt);
[StressN, Rn, lc_t] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
stress_0 = StressN + TV;


stress = stress_0;
devStress = Dev(stress);
tau = norm(devStress,'fro');
if tau > sigma0
    if eps0 == 0
        eps0 = strain_now;
        F_p_t_new = F_p_t;
    else
        gamDotP = a*abs(strain_now - eps0)^(b)*strainDot;
        F_p_t_new = (gamDotP/tau * (inv(F_ve_t)*devStress*F))*dt+F_p_t_dt;
    end
else
    eps0 = 0;
    F_p_t_new = F_p_t;
end

F_ve_t_new =  F * (F_p_t_new)^(-1);


F_e_t_new = F_ve_t_new * (F_b_t_new)^(-1);