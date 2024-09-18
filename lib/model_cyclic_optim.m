function [stress,trueStrain,exp_stress_updated,exp_strain_updated] = model_cyclic_optim(params,de,d0,dt,on,...
    expstrain,expstress,vp,params_F)
format longG 
F_i_t0 = eye(3); % Set F_i to 1 first
F_p_t0 = eye(3);
lc_max = 1;
d_t    = d0;
T      = zeros(3);
mu1 = params(1); mu2 = params(2); nu1 = params(3); nu2 = 0;
m = params(4); gamma_dot_0 = params(5); dG = params(6);
Ad = params(7); tau0 = params(8); d0s = params(9); m_tau = params(10);
a = params(11); b = params(12); sigma0 = params(13);
% kappa1 = 500;
% kappa2 = 500;
%kappa1 = 6 * mu1;
%kappa2 = 6 * mu2;
kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
kappa2 = 0;
params_A = [mu1, kappa1,nu1];
params_B = [mu2, kappa2,nu2];
X = 1 + 5*vp + 18*vp^2;

strain = 0;
F22    = 1; % Set F22 to 1 
TB_t = zeros(3);
TA_t = zeros(3);

strain_cycles = [0.016756727990000 0.019231870570000 0.021942198690000 ...
    0.025011823690000 0.028696511190000 0.033491720180000 0.035436974080000];
cycles = 6;
cycle = 1;
loading = 1;
unloading = 0;
timeVec = [0, dt];
trueStrain = [0, de];
stress = [0];
run = 1
params;
while run == 1
    if loading == 1
        time0 = timeVec(end-1);
        time1 = timeVec(end);
        F11 = 1 + trueStrain(end);
        
        ep  = 1e-4;
        iter_stress = 0; 
            while 1

                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;
                
                F = [F11 0 0; 0 F22 0; 0 0 F22]; % Incremental F at t+1
                J  = det(F);
                Fb = J^(-1/3) * F;
                F_b_t    = F_i_t0;
                F_b_t_dt = F_b_t; 
                F_p_t    = F_p_t0;
                F_p_t_dt = F_p_t;

                F_ve_t = Fb * (F_p_t)^(-1);
                F_e_t = F_ve_t*(F_b_t)^(-1);

                [TA_t, Rn, lc_t] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
                %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F


                %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F

                it_step = 0;
                while 1
                    % Do this loop until F_i_t_new is obtained

                    it_step = it_step + 1;


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F); % Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(F,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)

                    err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t]);
                    %err   = norm(F_b_t_new - F_b_t);
                    F_b_t = F_b_t_new;
                    F_e_t = F_e_t_new;
                    F_ve_t = F_ve_t_new;
                    F_p_t = F_p_t_new;
                    strain = strain_now;

                    if (err<1e-6)
                        break;
                    elseif it_step > 1000000
                        break;
                    else
                        break_pt = 1;
                    end
                end

                % After obtaining F_i_t_new we calculate the stress
                T = (1 - d_t) * (TA_t + TB_t);
                if isnan(F_b_t_new)
                    break;
                end

                % Get T_22 and if T_22 = 0 then start again with new F_i_t
                if (abs(T(2,2)) <= 0.1)
                    F_i_t0 = F_b_t;
                    F_p_t0 = F_p_t;
                    break;
                end

                Tnn = T(2,2);
                if Tn*Tnn < 0
                    ep = -ep/2;
                else
                    if abs(Tnn)>abs(Tn)
                        ep = -ep/2;
                    end
                end
                F22          = F22 - ep/2;
            end
         if (lc_t > lc_max)
             lc_max = lc_t;
             d_t    = 1 - exp(Ad * (1 - lc_max) );        
         end
         stress = [stress T(1,1)];
         if trueStrain(end) > strain_cycles(cycle)
             loading = 0;
             unloading = 1;
             trueStrain = [trueStrain trueStrain(end)-de];
             timeVec = [timeVec timeVec(end)+dt];
             if cycle == cycles
                run = 0;
                trueStrain(end) = [];
             end
         else
             loading = 1;
             unloading = 0;
             trueStrain = [trueStrain trueStrain(end)+de];
             timeVec = [timeVec timeVec(end)+dt];
         end
    else
        time0 = timeVec(end-1);
        time1 = timeVec(end);
        F11 = 1 + trueStrain(end);
        ep  = 1e-4;
        iter_stress = 0; 
             while 1

                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;

                F = [F11 0 0; 0 F22 0; 0 0 F22]; % Incremental F at t+1
                J  = det(F);
                Fb = J^(-1/3) * F;                
                F_b_t    = F_i_t0;
                F_b_t_dt = F_b_t; 
                F_p_t    = F_p_t0;
                F_p_t_dt = F_p_t;

                F_ve_t = Fb * (F_p_t)^(-1);
                F_e_t = F_ve_t*(F_b_t)^(-1);

                [TA_t, Rn, lc_t] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
                %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F


                %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F

                it_step = 0;
                while 1
                    % Do this loop until F_i_t_new is obtained

                    it_step = it_step + 1;


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F);% Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(F,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)

                    err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t]);
                    %err   = norm(F_b_t_new - F_b_t);
                    F_b_t = F_b_t_new;
                    F_e_t = F_e_t_new;
                    F_ve_t = F_ve_t_new;
                    F_p_t = F_p_t_new;
                    strain = strain_now;

                    if (err<1e-6)
                        break;
                    elseif it_step > 1000000
                        break;
                    else
                        break_pt = 1;
                    end
                end

                % After obtaining F_i_t_new we calculate the stress
                T = (1 - d_t) * (TA_t + TB_t);
                if isnan(F_b_t_new)
                    break;
                end

                % Get T_22 and if T_22 = 0 then start again with new F_i_t
                if (abs(T(2,2)) <= 0.1)
                    F_i_t0 = F_b_t;
                    F_p_t0 = F_p_t;
                    break;
                end

                Tnn = T(2,2);
                if Tn*Tnn < 0
                    ep = -ep/2;
                else
                    if abs(Tnn)>abs(Tn)
                        ep = -ep/2;
                    end
                end
                F22          = F22 + ep/2;
            end
         if (lc_t > lc_max)
             lc_max = lc_t;
             d_t    = 1 - exp(Ad * (1 - lc_max) );        
         end
         stress = [stress T(1,1)];
         if stress(end) < 1e-4
             loading = 1;
             unloading = 0;
             trueStrain = [trueStrain trueStrain(end)+de];
             timeVec = [timeVec timeVec(end)+dt];
             cycle = cycle + 1;
         else
             loading = 0;
             unloading = 1;
             trueStrain = [trueStrain trueStrain(end)-de];
             timeVec = [timeVec timeVec(end)+dt];
         end
    end                
end
plot(trueStrain,stress);
pause(3);

close all;

end