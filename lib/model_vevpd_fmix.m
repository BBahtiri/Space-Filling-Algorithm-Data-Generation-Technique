function [stress,F_tot,F_ve_tot,F_v_tot,F_e_tot,F_vp_tot,J_tot,timeVec,...
    T_model_all,T_model_networkA,T_model_networkB,T_model_volumetric,d_total,...
    strain_gl,tau_tot,trueStrain1,b_alm,C_all,rdn] = model_calibr(params,de,d0,dt,vp,params_F,nnetwork,...
    cyclic_loading,et_load,et_unload,trueStrain1,trueStrain2,trueStrain3,trueStrain12,...
    trueStrain13,trueStrain21,trueStrain23,trueStrain31,trueStrain32,perturbation)
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
alphaZ = params_F(5);
alphaT = params_F(3);
X = (1 + 5*vp + 18*vp^2)*alphaZ*alphaT;
% kappa1 = 500;
% kappa2 = 500;
%kappa1 = 6 * mu1;
%kappa2 = 6 * mu2;
kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
kappa2 = 0;
params_A = [mu1, kappa1,nu1];
params_B = [mu2, kappa2,nu2];
F22    = 1; % Set F22 to 1 
T_model_loading = zeros(length(trueStrain1),1);
T_model_networkA = {};
T_model_networkB = {};
T_model_volumetric = {};
T_model_all = {};
strain_gl = {};
b_alm = {};
b_hist = {};
d_total         = zeros(length(trueStrain1),1);
J_tot         = zeros(length(trueStrain1),1);
t_tot = zeros(length(trueStrain1),1);
F_tot = {};
F_ve_t = {};
F_v_t = {};
F_e_t = {};
F_vp_t = {};
F_vp_tot = {};
F_ve_tot = {};
F_v_tot = {};
F_e_tot = {};
timeVec = zeros(length(trueStrain1),1);
tau_tot = zeros(length(trueStrain1),1);
gam_tot = zeros(length(trueStrain1),1);
stress = [];
TB_t = zeros(3);
strain = 0;
C_all = {};
C_num= zeros(9);
epsilon = 1e-4;
rdn = randi(9,1);
% for i = 1:9
%     switch i
%         case 1
% 	    if i == rdn
%             	trueStrain1(end+1) = trueStrain1(end)+epsilon;
% 	    else
% 		trueStrain1(end+1) = trueStrain1(end);
% 	    end
%         case 2
% 	    if i == rdn
%             	trueStrain2(end+1) = trueStrain2(end)+epsilon;
% 	    else
% 		trueStrain2(end+1) = trueStrain2(end);
% 	    end
%         case 3 
% 	    if i == rdn
%             	trueStrain3(end+1) = trueStrain3(end)+epsilon;
% 	    else
% 		trueStrain3(end+1) = trueStrain3(end);
% 	    end
%         case 4
% 	    if i == rdn
%             	trueStrain12(end+1) = trueStrain12(end)+epsilon;
% 	    else
% 		trueStrain12(end+1) = trueStrain12(end);
% 	    end
%         case 5
% 	    if i == rdn
%             	trueStrain13(end+1) = trueStrain13(end)+epsilon;
% 	    else
% 		trueStrain13(end+1) = trueStrain13(end);
% 	    end
%         case 6
% 	    if i == rdn
%             	trueStrain21(end+1) = trueStrain21(end)+epsilon;
% 	    else
% 		trueStrain21(end+1) = trueStrain21(end);
% 	    end
%         case 7
% 	    if i == rdn
%             	trueStrain23(end+1) = trueStrain23(end)+epsilon;
% 	    else
% 		trueStrain23(end+1) = trueStrain23(end);
% 	    end
%         case 8
% 	    if i == rdn
%             	trueStrain31(end+1) = trueStrain31(end)+epsilon;
% 	    else
% 		trueStrain31(end+1) = trueStrain31(end);
% 	    end
%         case 9
% 	    if i == rdn
%             	trueStrain32(end+1) = trueStrain32(end)+epsilon;
% 	    else
% 		trueStrain32(end+1) = trueStrain32(end);
% 	    end
%     end  
% end

if cyclic_loading == 0
for j = 1:length(trueStrain1)
        F11 = trueStrain1(j); % Increment for F11
        F22 = trueStrain2(j);
        F33 = trueStrain3(j);
        F12 = trueStrain12(j);
        F13 = trueStrain13(j);
        F21 = trueStrain21(j);
        F23 = trueStrain23(j);
        F31 = trueStrain31(j);
        F32 = trueStrain32(j);
        F = [F11 F12 F13; F21 F22 F23; F31 F32 F33]; % Incremental F at t+1
        J  = det(F);
        Fb = J^(-1/3) * F;
        Cb = F'*F;
        bF = F*F';
        
        E_gl = 1/2*(Cb-eye(3));
        F_b_t    = F_i_t0;
        F_b_t_dt = F_b_t; 
        F_p_t    = F_p_t0;
        F_p_t_dt = F_p_t;
        
        F_ve_t = Fb * (F_p_t)^(-1);
        F_e_t = F_ve_t*(F_b_t)^(-1);
        
        [TA_t, Rn, lc_t,T_vol] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
            it_step = 0;
           
              while 1
                % Do this loop until F_i_t_new is obtained

                it_step = it_step + 1;


                [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F); % Stress for elastic spring with F_e_t coming from F_i_t

                [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,strain_now] = Fi_backward(Fb,F_b_t,F_b_t_dt,...
                    F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)     
                err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t],'fro');
                F_b_t = F_b_t_new;
                F_e_t = F_e_t_new;
                F_ve_t = F_ve_t_new;
                F_p_t = F_p_t_new;


                if(it_step > 1e3)
                    break;
                end
                if (err<1e-6)
                    break;
                else
                    break_pt = 1;
                end
              end
              it_step_all{j} = it_step;
            
            % After obtaining F_i_t_new we calculate the stress
            Td = TA_t + TB_t;
            T = (1 - d_t) * (TA_t + TB_t);

            F_i_t0 = F_b_t;
            F_p_t0 = F_p_t;
            strain = strain_now;
            if isnan(F_b_t)
                break;
            end

            if perturbation == 1
                taustress = T.*J;
                C_num= zeros(9);
                alpha = 1e-4;
                for ip = 1:3
                    
                    ini = 3*(ip-1)+1:3*ip;
                    ei = zeros(1,3);
                    ei(ip) = 1;
                    
                    for jp = 1:3
                        
                        ej = zeros(1,3);
                        ej(jp) = 1;
                        
                        % perturbed deformation gradient (Fb_ij)
                        delF = (alpha/2)*(ei'*ej*F + ej'*ei*F);
                        Fp   = F + delF;
                        
                        % perturbed second Piola-Kirchhoff stress
                        invcp      = inv(Fp'*Fp);
                        Jp   = det(Fp);
                        lnJp = log(Jp);
                        be= Fp*Fp';
                        [cauchyp] = rheo_model(Fp,Jp,be,params,params_F,de,dt,F_i_t0,F_p_t0,strain,...
                            d_t,eps0);
                        taup = cauchyp*Jp;
                               
                        % calculate dtau/dD
                        inj = 3*(jp-1)+1:3*jp;        
                        C_num(ini,inj) = (taup - taustress) / alpha;
                        
                    end
                end
            end
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
    b_alm{j,1} = bF;
    timeVec(j+1) = timeVec(j)+dt;
    F_tot{j,1 } = F;
    F_ve_tot{j,1} = F_ve_t;
    F_v_tot{j,1} = F_b_t;
    F_e_tot{j,1} = F_e_t;
    F_vp_tot{j,1} = F_p_t;
    J_tot(j) = J;
    stress = T_model_loading;
    C_all{j} = C_num;
end
elseif cyclic_loading == 1
    strain_cycles = [0.0087 0.012 0.0155 ...
    0.0192 0.022 0.027 0.032];
% % 
%     strain_cycles = [0.0025 0.005 0.0075 ...
%      0.01 0.0125 0.015 0.0175];
%     strain_cycles = [0.005 0.01 0.015 ...
%     0.02 0.025 0.03 0.035];
%     strain_cycles = [0.01 0.02 0.03 ...
%     0.04 0.05 0.06 0.07];
    cycles = 7;
    cycle = 1;
    loading = 1;
    unloading = 0;
    timeVec = [0, dt];
%     dt = 0.1;
    de = single(5e-5);
    de2 = single(6e-5);
    de3 = single(3e-5);
    dt = de/et_load;
    trueStrain = [0, de];
    trueStrain2 = [0, de2];
    trueStrain3 = [0, de3];
    stress = [0];
    run = 1;
    j = 1;
    
while run == 1
    if loading == 1
        timeVec = [0, dt];
        time0 = timeVec(end-1);
        time1 = timeVec(end);
        F11 = 1 + trueStrain(end);
%         F22 = 1.0 + trueStrain2(end);
%         F33 = 1.0 + trueStrain3(end);
        F22 = 1.0;
        F33 = 1.0;
        ep  = 1e-4;
        iter_stress = 0; 

                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;
                
                F = [F11 0 0; 0 F22 0; 0 0 F33]; % Incremental F at t+1
                J  = det(F);
                Fb = J^(-1/3) * F;
                Cb = F'*F;
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
                T = (1 - d_t) * (TA_t + TB_t);
                Td = TA_t + TB_t;
                if isnan(F_b_t_new)
                    break;
                end
                F_i_t0 = F_b_t;
                F_p_t0 = F_p_t;
                strain = strain_now;
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
            F_tot{j,1 } = F;
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
             de = single(5e-5);
             de2 = single(6e-5);
             de3 = single(3e-5);
             dt = de/et_load;
             trueStrain = [trueStrain trueStrain(end)-de];
             trueStrain2 = [trueStrain2 trueStrain2(end)-de2];
             trueStrain3 = [trueStrain3 trueStrain3(end)-de3];
             timeVec = [timeVec timeVec(end)+dt];
             if cycle == cycles
                run = 0;
                trueStrain(end) = [];
             end
         else
             loading = 1;
             unloading = 0;
             trueStrain = [trueStrain trueStrain(end)+de];
             trueStrain2 = [trueStrain2 trueStrain2(end)+de2];
             trueStrain3 = [trueStrain3 trueStrain3(end)+de3];
             timeVec = [timeVec timeVec(end)+dt];
         end
    else
        time0 = timeVec(end-1);
        time1 = timeVec(end);
        F11 = 1 + trueStrain(end);
%         F22 = 1.0 + trueStrain2(end);
%         F33 = 1.0 + trueStrain3(end);
        F22 = 1.0;
        F33 = 1.0;
        ep  = 1e-4;
        iter_stress = 0; 
                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;

                F = [F11 0 0; 0 F22 0; 0 0 F33]; % Incremental F at t+1
                J  = det(F);
                Fb = J^(-1/3) * F;  
                Cb = F'*F;
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
                T = (1 - d_t) * (TA_t + TB_t);
                Td = TA_t + TB_t;
                if isnan(F_b_t_new)
                    break;
                end
                F_i_t0 = F_b_t;
                F_p_t0 = F_p_t;
                strain = strain_now;

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
            F_tot{j,1 } = F;
            F_ve_tot{j,1} = F_ve_t;
            F_v_tot{j,1} = F_b_t;
            F_e_tot{j,1} = F_e_t;
            F_vp_tot{j,1} = F_p_t;
            J_tot(j) = J;
            j = j + 1;
         if stress(end) < 1e-4
             loading = 1;
             unloading = 0;
             de = single(5e-5);
             de2 = single(6e-5);
             de3 = single(3e-5);
             dt = de/et_load;
             trueStrain = [trueStrain trueStrain(end)+de];
             trueStrain2 = [trueStrain2 trueStrain2(end)+de2];
             trueStrain3 = [trueStrain3 trueStrain3(end)+de3];
             timeVec = [timeVec timeVec(end)+dt];
             cycle = cycle + 1;
         else
             loading = 0;
             unloading = 1;
             trueStrain = [trueStrain trueStrain(end)-de];
             trueStrain2 = [trueStrain2 trueStrain2(end)-de2];
             trueStrain3 = [trueStrain3 trueStrain3(end)-de3];
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

% Function to calculate stress by using rheological model
function [cauchy] = rheo_model(F,J,b,params,params_F,de,dt,...
    F_i_t0,F_p_t0,strain,d_t,eps0)

    % wnp/wgf = 0/0 - 0.1/0 - 0.2/0 - 0/0.6 - 0.15/0.45 - 0.1/0.3 - 0.15/0.15
    wnp   = 0.0;              % nanoparticle weight fraction
    wgf   = 0.0;              % glass fiber weight fraction
    wp    = 1 - wnp - wgf;    % polymer weight fraction

    % densities
    ro_p  = 1.20;               % density of polymer (g/ml)
    ro_np = 3.00;               % density of nanoparticle (g/ml)
    ro_gf = 2.55;               % density of glass fiber (g/ml)

    % nanoparticle volume fraction
    vp    = wnp * ro_p / (ro_np + wnp*ro_p - ro_np*wnp);
%     eps0 = 0;
%     F_i_t0 = eye(3); % Set F_i to 1 first
%     F_p_t0 = eye(3);
    cauchy = zeros(3);
    mu1 = params(1); mu2 = params(2); nu1 = params(3); nu2 = 0;
    m = params(4); gamma_dot_0 = params(5); dG = params(6);
    Ad = params(7); tau0 = params(8); d0s = params(9); m_tau = params(10);
    a = params(11); b = params(12); sigma0 = params(13);
    y0 = params(14); x0 = params(15); a_t = params(16); b_t = params(17);
    X = 1 + 5*vp + 18*vp^2;
    kappa1 = (2*mu1*(1+nu1))/(3*(1-2*nu1));
    kappa2 = 0;
    params_A = [mu1, kappa1,nu1];
    params_B = [mu2, kappa2,nu2];
    Fb = J^(-1/3) * F;
    
    
    F_b_t    = F_i_t0;
    F_b_t_dt = F_b_t; 
    F_p_t    = F_p_t0;
    F_p_t_dt = F_p_t;
    
    F_ve_t = Fb * (F_p_t)^(-1);
    F_e_t = F_ve_t*(F_b_t)^(-1);
    [TA_t, Rn, lc_t,T_vol] = NH_3D(F_ve_t,params_A,1,J,X,params_F);
    it_step = 0;
%     strain = 0;
     while 1
        % Do this loop until F_i_t_new is obtained

        it_step = it_step + 1;


        [TB_t,R,lc_B]              = NH_3D(F_e_t,params_B,0,J,X,params_F); % Stress for elastic spring with F_e_t coming from F_i_t

        [F_e_t_new,F_b_t_new,F_p_t_new,F_ve_t_new,tau,gamma,~] = Fi_backward(Fb,F_b_t,F_b_t_dt,...
            F_p_t,F_p_t_dt,F_ve_t,TB_t,R,params,d_t,TA_t,dt,de,J,X,params_F,strain); % Obtaining F_e_t_new from dashpot (time integration)     
        err   = norm([F_b_t_new; F_p_t_new] - [F_b_t; F_p_t],'fro');
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
    
    cauchy = (1 - d_t) * (TA_t + TB_t);
     
end

end

