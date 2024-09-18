function [stress,trueStrain] = model_cyclic(params,de,d0,dt,on)
format long 
F_i_t0 = eye(3); % Set F_i to 1 first
F_p_t0 = eye(3);
lc_max = 1;
d_t    = d0;
T      = zeros(3);
mu1 = params(1); mu2 = params(2); nu1 = params(3); nu2 = params(4);
m = params(5); gamma_dot_0 = params(6); dG = params(7);
Ad = params(8); tau0 = params(9); d0s = params(10); m_tau = params(11);
a = params(12); b = params(13); sigma0 = params(14);

% kappa1 = 500;
% kappa2 = 500;
%kappa1 = 6 * mu1;
%kappa2 = 6 * mu2;
kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
kappa2 = (2*mu2*(1+nu2))/(3*(1-2*nu2));
params_A = [mu1, kappa1,nu1];
params_B = [mu2, kappa2,nu2];
F22    = 1; % Set F22 to 1 
TB_t = zeros(3);
TA_t = zeros(3);

strain_cycles = [0.016756727990000 0.019231870570000 0.021942198690000 ...
    0.025011823690000 0.028696511190000 0.033491720180000 0.035436974080000];
cycles = length(strain_cycles);
cycle = 1;
loading = 1;
unloading = 0;
timeVec = [0, dt];
trueStrain = [0, de];
stress = [0];
run = 1;
params
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


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J); % Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma] = Fi_backward(F,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J); % Obtaining F_e_t_new from dashpot (time integration)

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
                run = 0
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

                [TA_t, Rn, lc_t] = NH_3D(Fb,params_A,1,J);
                %[TA_t] = stress_N(F,params_TA,vf,a0,vf_ratio,vnpp); % Stress for hyperealstic spring with F


                %F_e_t    = F * (F_b_t)^(-1); % F_e_t from F_i_t and F

                it_step = 0;
                while 1
                    % Do this loop until F_i_t_new is obtained

                    it_step = it_step + 1;


                    [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J);% Stress for elastic spring with F_e_t coming from F_i_t

                    [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma] = Fi_backward(F,F_b_t,F_b_t_dt,...
                        F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de); % Obtaining F_e_t_new from dashpot (time integration)

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

if on == 1
    [pks_max_calc,locs_max_calc] = findpeaks(trueStrain);
    [pks_min_calc,locs_min_calc] = findpeaks(-trueStrain);
    locs_max_calc = locs_max_calc';
    locs_min_calc = locs_min_calc';
    blocks_calc = sort([1; locs_max_calc; locs_min_calc]);
    interp_stress_calc = [];
    new_strain = [];
    k = 1;
    for i = 2:length(blocks_calc)
        block_curr_str = trueStrain(blocks_calc(i-1):blocks_calc(i));
        block_curr_stress = stress(blocks_calc(i-1):blocks_calc(i));
        block_curr_calc_interp = interp1(block_curr_str, block_curr_stress, block_curr_str);
        interp_stress_calc = [interp_stress_calc block_curr_calc_interp];
        new_strain = [new_strain block_curr_str];
        k = k + 1;

    end
    % Last loading
    block_curr_str = trueStrain(blocks_calc(end):end);
    block_curr_stress = stress(blocks_calc(end):end);
    block_curr_calc_interp = interp1(block_curr_str, block_curr_stress, block_curr_str);
    interp_stress_calc = [interp_stress_calc block_curr_calc_interp];
    new_strain = [new_strain block_curr_str];
    stress = interp_stress_calc;
    trueStrain = new_strain;
    plot(trueStrain,stress);
    pause(3);
end

end