function [stress,F_tot,F_ve_tot,F_v_tot,F_e_tot,F_vp_tot,J_tot,timeVec,...
    T_model_all,T_model_networkA,T_model_networkB,T_model_volumetric,d_total,...
    strain_gl,tau_tot] = model_calibr(params,de,d0,dt,trueStrain,vp,params_F,nnetwork,...
    cyclic_loading,et_load,et_unload)
format longG
global eps0
eps0 = 0;
F_i_t0 = eye(3); % Set F_i to 1 first
F_p_t0 = eye(3);
lc_max = 1;
d_t    = d0;
T      = zeros(3);
mu1 = params(1); mu2 = params(2); nu1 = params(3); nu2 = 0;
m = params(4); gamma_dot_0 = params(5); dG = params(6);
Ad = params(7); tau0 = params(8); d0s = params(9); m_tau = params(10);
a = params(11); b = params(12); sigma0 = params(13);
y0 = params(14); x0 = params(15); a_t = params(16); b_t = params(17);
zita   = params_F(4);
alphaZ = params_F(5);
X = (1 + 5*vp + 18*vp^2);
% kappa1 = 500;
% kappa2 = 500;
%kappa1 = 6 * mu1;
%kappa2 = 6 * mu2;
kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
kappa2 = 0;
params_A = [mu1, kappa1,nu1];
params_B = [mu2, kappa2,nu2];
F22    = 1; % Set F22 to 1 
T_model_loading = zeros(length(trueStrain),1);
T_model_networkA = {};
T_model_networkB = {};
T_model_volumetric = {};
T_model_all = {};
strain_gl = {};
d_total         = zeros(length(trueStrain),1);
J_tot         = zeros(length(trueStrain),1);
t_tot = zeros(length(trueStrain),1);
F_tot = {};
F_ve_t = {};
F_v_t = {};
F_e_t = {};
F_vp_t = {};
timeVec = zeros(length(trueStrain),1);
tau_tot = zeros(length(trueStrain),1);
gam_tot = zeros(length(trueStrain),1);

TB_t = zeros(3);
strain = 0;

if cyclic_loading == 0
for j = 1:length(trueStrain)
    F11 = 1 + trueStrain(j); % Increment for F11
    ep  = 1e-5;

    % loop to impose conditions of pure tensile loading (i.e., T_yy = T_zz = 0)
    iter_stress = 0;
    while 1
        
        Tn = T(2,2); % Stress in uniaxial direction
        iter_stress = iter_stress + 1;
        
        F = [F11 0 0; 0 F22 0; 0 0 F22]; % Incremental F at t+1
        J  = det(F);
        Fb = J^(-1/3) * F;
        Cb = Fb'*Fb;
        E_gl = 1/2*(Cb-eye(3));
        F_b_t    = F_i_t0;
        F_b_t_dt = F_b_t; 
        F_p_t    = F_p_t0;
        F_p_t_dt = F_p_t;
        
        F_ve_t = Fb * (F_p_t)^(-1);
        F_e_t = F_ve_t*(F_b_t)^(-1);
        
        [TA_t, Rn, lc_t,T_vol] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
        %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F
        
        
        %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F
        
            it_step = 0;
              while 1
                % Do this loop until F_i_t_new is obtained

                it_step = it_step + 1;


                [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F); % Stress for elastic spring with F_e_t coming from F_i_t

                [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(Fb,F_b_t,F_b_t_dt,...
                    F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)     
                err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t],'fro');
                %err   = norm(F_b_t_new - F_b_t);
                F_b_t = F_b_t_new;
                F_e_t = F_e_t_new;
                F_ve_t = F_ve_t_new;
                F_p_t = F_p_t_new;


                if(it_step > 1e5)
                    it_step
                end
                if (err<1e-6)
                    break;
                else
                    break_pt = 1;
                end
            end

            % After obtaining F_i_t_new we calculate the stress
            Td = TA_t + TB_t;
            T = (1 - d_t) * (TA_t + TB_t);
            if isnan(F_b_t_new)
                break;
            elseif it_step > 10000
                break;
            end

            % Get T_22 and if T_22 = 0 then start again with new F_i_t
            if (abs(T(2,2)) <= 0.1)
                F_i_t0 = F_b_t;
                F_p_t0 = F_p_t;
                strain = strain_now;
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
            if isnan(F_b_t_new)
               break;
            elseif it_step > 10000
               break;
            end
        %     
            % update damage
             if (lc_t > lc_max)
                 lc_max = lc_t;
                d_t    = 1 - exp(Ad * (1 - lc_max) );        
             end

            Tdev                 = Dev(T);
            %T_model_loading(j+1) = sqrt(3/2 *trace(Tdev*Tdev'));          % Von Mises Stress
            d_total(j)        = d_t;
            tau_tot(j)        = tau;
            gam_tot(j)        = gamma;
            T_model_loading(j) = T(1,1);

    T_model_all{j,1} = Td;
    T_model_networkA{j,1} = TA_t;
    T_model_networkB{j,1} = TB_t;
    T_model_volumetric{j,1} = T_vol;
    strain_gl{j,1} = E_gl;
    timeVec(j+1) = timeVec(j)+dt;
    F_tot{j,1 } = Fb;
    F_ve_tot{j,1} = F_ve_t;
    F_v_tot{j,1} = F_b_t;
    F_e_tot{j,1} = F_e_t;
    F_vp_tot{j,1} = F_p_t;
    J_tot(j) = J;
    stress = T_model_loading;
end
elseif cyclic_loading == 1
    strain_cycles = [0.0087 0.012 0.0155 ...
    0.0192 0.022 0.027 0.032];
    cycles = 7;
    cycle = 1;
    loading = 1;
    unloading = 0;
    timeVec = [0, dt];
    de = et_load * dt;
    trueStrain = [0, de];
    stress = [0];
    run = 1;
    j = 1;
    
while run == 1
    if loading == 1
        timeVec = [0, dt];
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
                Cb = Fb'*Fb;
                E_gl = 1/2*(Cb-eye(3));
                F_b_t    = F_i_t0;
                F_b_t_dt = F_b_t; 
                F_p_t    = F_p_t0;
                F_p_t_dt = F_p_t;

                F_ve_t = Fb * (F_p_t)^(-1);
                F_e_t = F_ve_t*(F_b_t)^(-1);

                [TA_t, Rn, lc_t,T_vol] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
                %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F


                %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F

                it_step = 0;
                while 1
                    % Do this loop until F_i_t_new is obtained

                    it_step = it_step + 1;


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F); % Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(Fb,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain);% Obtaining F_e_t_new from dashpot (time integration)

                    err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t]);
                    %err   = norm(F_b_t_new - F_b_t);
                    F_b_t = F_b_t_new;
                    F_e_t = F_e_t_new;
                    F_ve_t = F_ve_t_new;
                    F_p_t = F_p_t_new;

                    if (err<1e-6)
                        break;
                    elseif it_step > 1000000
                        break;
                    else
                        break_pt = 1;
                    end
                end

                % After obtaining F_i_t_new we calculate the stress
                Td = TA_t + TB_t;
                T = (1 - d_t) * (TA_t + TB_t);
                if isnan(F_b_t_new)
                    break;
                end

                % Get T_22 and if T_22 = 0 then start again with new F_i_t
                if (abs(T(2,2)) <= 0.1)
                    F_i_t0 = F_b_t;
                    F_p_t0 = F_p_t;
                    strain = strain_now;
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
            T_model_all{j,1} = Td;
            T_model_networkA{j,1} = TA_t;
            T_model_networkB{j,1} = TB_t;
            T_model_volumetric{j,1} = T_vol;
            strain_gl{j,1} = E_gl;
            %timeVec(j+1) = timeVec(j)+dt;
            F_tot{j,1 } = Fb;
            F_ve_tot{j,1} = F_ve_t;
            F_v_tot{j,1} = F_b_t;
            F_e_tot{j,1} = F_e_t;
            F_vp_tot{j,1} = F_p_t;
            J_tot(j) = J;
            j = j + 1;
         if trueStrain(end) > strain_cycles(cycle)
             loading = 0;
             unloading = 1;
             timeVec = [0, dt];
             
             de = et_unload * dt;
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
                Cb = Fb'*Fb;
                E_gl = 1/2*(Cb-eye(3));                
                F_b_t    = F_i_t0;
                F_b_t_dt = F_b_t; 
                F_p_t    = F_p_t0;
                F_p_t_dt = F_p_t;

                F_ve_t = Fb * (F_p_t)^(-1);
                F_e_t = F_ve_t*(F_b_t)^(-1);

                [TA_t, Rn, lc_t,T_vol] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
                %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F


                %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F

                it_step = 0;
                while 1
                    % Do this loop until F_i_t_new is obtained

                    it_step = it_step + 1;


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F);% Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(Fb,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)

                    err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t]);
                    %err   = norm(F_b_t_new - F_b_t);
                    F_b_t = F_b_t_new;
                    F_e_t = F_e_t_new;
                    F_ve_t = F_ve_t_new;
                    F_p_t = F_p_t_new;
                    

                    if (err<1e-6)
                        break;
                    elseif it_step > 1000000
                        break;
                    else
                        break_pt = 1;
                    end
                end

                % After obtaining F_i_t_new we calculate the stress
                Td = TA_t + TB_t;
                T = (1 - d_t) * (TA_t + TB_t);
                if isnan(F_b_t_new)
                    break;
                end

                % Get T_22 and if T_22 = 0 then start again with new F_i_t
                if (abs(T(2,2)) <= 0.1)
                    F_i_t0 = F_b_t;
                    F_p_t0 = F_p_t;
                    strain = strain_now;
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
            T_model_all{j,1} = Td;
            T_model_networkA{j,1} = TA_t;
            T_model_networkB{j,1} = TB_t;
            T_model_volumetric{j,1} = T_vol;
            strain_gl{j,1} = E_gl;
            %timeVec(j+1) = timeVec(j)+dt;
            F_tot{j,1 } = Fb;
            F_ve_tot{j,1} = F_ve_t;
            F_v_tot{j,1} = F_b_t;
            F_e_tot{j,1} = F_e_t;
            F_vp_tot{j,1} = F_p_t;
            J_tot(j) = J;
            j = j + 1;
         if stress(end) < 1e-4
             loading = 1;
             unloading = 0;
             de = et_load * dt;
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
end



function [y] = wrapper(matrix)
    F11 = matrix(1,1);
    F12 = matrix(1,2);
    F13 = matrix(1,3);
    F21 = matrix(2,1);
    F22 = matrix(2,2);
    F23 = matrix(2,3);
    F31 = matrix(3,1);
    F32 = matrix(3,2);
    F33 = matrix(3,3);
    y = [F11 F12 F13 F21 F22 F23 F31 F32 F33];
end

end

