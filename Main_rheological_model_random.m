clear all;
clc;
close all;
global a1 a2 a3

p_s = [pwd '/lib']
addpath(p_s)
addpath('C:/Users/bahtiri/Documents/postprocessingMD/Matlab/matlab2tikz/matlab2tikz-master/src/')


% wnp/wgf = 0/0 - 0.1/0 - 0.2/0 - 0/0.6 - 0.15/0.45 - 0.1/0.3 - 0.15/0.15
wnp   = 0.0;              % nanoparticle weight fraction
wgf   = 0.0;              % glass fiber weight fraction
wp    = 1 - wnp - wgf;    % polymer weight fraction


% densities
ro_p  = 1.20;               % density of polymer (g/ml)
ro_np = 3.00;               % density of nanoparticle (g/ml)
ro_gf = 2.55;               % density of glass fiber (g/ml)

% nanoparticle volume fraction
vnp    = wnp * ro_p / (ro_np + wnp*ro_p - ro_np*wnp);

% cyclic_exp_full = 0;
% cyclic = 0;
% if cyclic_exp_full == 1
% %     filename = "exp_data_cyclic/epoxy_1.txt";
% %     [timeVec, trueStrain, expStress] = read_in(filename,1); 
% %     p = polyfit(timeVec(100:200),trueStrain(100:200),1);
%     et_load = 0.00015;
%     et_unload = 0.00030;
%     dt = 0.1;
%     de = et * dt;
% else
%     filename = ["monotonic_edited.txt"] ;
%     [timeVec, trueStrain_mono, expStress] = read_in(filename,0); 
%     p = polyfit(timeVec(100:200),trueStrain_mono(100:200),1);
% %     et = 5e-4; % et3
% %     et = 1e-3; % et2
%     et = 9e-4; % et1
% %   et = 1.6e-3;
% %    dt_mono = 0.1;
%     ef = 0.06;
%     de_mono = mean(abs(diff(trueStrain_mono)));
%     de_mono = 5e-5;
%     dt_mono = de_mono/et;
%     strain_step  = round(ef/de_mono); 
%     trueStrain_mono = linspace(0, ef, strain_step);
% %     et_load =5e-4; %et3
% %     et_unload = 5e-4; %et3
% %     et_load =1e-3; %et2
% %     et_unload = 1e-3; %et2
%     et_load =9e-4; %et1
%     et_unload =9e-4; %et1
% %     et_load = 0.00015;
% %     et_unload = 0.0003;
% end

%TEST
mu1 = 760;
mu2 = 790;
nu1 = 0.23;
nu2 = 0.0;
m = 0.657;
gamma_dot_0 = 9.7746e11;
dG = 1.9761e-19;
Ad = 320;
tau0 = 85;
d0s = 0.1;
m_tau = 2;
y0 = 75;
x0 = 0.2369;
a_t = -48.23;
b_t = 0.06786;
a = 0.179;
% b = 0.910;
sigma0 =5.5;
Temper = 296;
% zita = 0.011;
zita = 0.0;
alphaZ = 1 + 0.057*zita.^2-9.5.*zita;
% alphaZ = 0.039;
Tref = 296;
alphaT = 2 - exp(0.01093*(Temper - Tref));
b = 0.910;
params = [mu1 mu2 nu1 m gamma_dot_0 dG Ad tau0 d0s m_tau a b sigma0 y0 x0 ...
    a_t b_t];
params_F = [Temper Tref alphaT zita alphaZ wnp];
nnetwork = 1;
d0 = 0;
de_mono = 5e-5;
et = 1e-3;
dt_mono = de_mono/et;
et_load =1e-3;
et_unload = 1e-3;


%% Random simulations
sim_strain = {};
peaks = 5;
% set(0,'DefaultFigureVisible','off');
target_fold = sprintf('simulations_const_random_%d',peaks);
mkdir(target_fold);
target_fold2 = sprintf('simulations_const_random_pert_%d',peaks);
mkdir(target_fold2);
p = haltonset(9);
p = scramble(p,'RR2');
X0 = p(1:1000,:);
M1 = normalize(X0(:,1),'range',[0.9,1.1]);
M2 = normalize(X0(:,2),'range',[0.9,1.1]);
M3 = normalize(X0(:,3),'range',[0.9,1.1]);
M4 = normalize(X0(:,4),'range',[-0.05,0.05]);
M5 = normalize(X0(:,5),'range',[-0.05,0.05]);
M6 = normalize(X0(:,6),'range',[-0.05,0.05]);
M7 = normalize(X0(:,7),'range',[-0.05,0.05]);
M8 = normalize(X0(:,8),'range',[-0.05,0.05]);
M9 = normalize(X0(:,9),'range',[-0.05,0.05]);
f1 = figure('Name','Uniform random');
scatter3(M1,M2,M3) 
hold on;
% hold on;
% plot([0.93 0.93],[0.93 1.07],'LineWidth',2,'Color','r');
% plot([1.07 1.07],[0.93 1.07],'LineWidth',2,'Color','r');
% plot([0.93 1.07],[0.93 0.93],'LineWidth',2,'Color','r');
% plot([0.93 1.07],[1.07 1.07],'LineWidth',2,'Color','r');
% plot([1.0 1.0],[0.93 1.07],'LineWidth',2,'Color','b');
% plot([0.93 1.07],[1.0 1.0],'LineWidth',3,'Color','b');
ylim([0.9 1.1])
xlim([0.9 1.1])
start = ones(length(X0),1);
origin = [1;1;1];
startin_ptx = 1.0;
startin_pty = 1.0;
startin_ptz = 1.0;
lpd = ones(3,peaks);
lps = zeros(6,peaks);
combination = 123;
f2 = figure('Name','Rheology model');
hold on;
f4 = figure('Name','Rheology shear stress');
f3 = figure('Name','Deformation Gradient');
validate = 0;
pos_comb = [123 1 2 3 12 13 23 9 112 113 121 1230 131 132 ...
            212 213 221 223 231 232 312 313 321 323 331 332 9 9 9 9 9 9 9];
%pos_comb = [123 1 2 3 12 13 23];
% pos_comb = [12];
% pos_comb = 9;
wnp_m = [0;0.1];
zitam = [0;0.01];
temps = [296;333];
time_cons = zeros(100,4);
time_ml = zeros(100,4);
if validate == 0
    for sims = 1
        ind = randperm(numel(pos_comb), 1);
        combination =  pos_comb(ind) %Randomly select combination of F
        %Starting points
        startin_ptx = 1.0; 
        startin_pty = 1.0;
        startin_ptz = 1.0;
        % Select randomly states for the simulation from the uniform points
        next = randi([1 length(X0)],peaks,1);
        statexx = M1(next);
        stateyy = M2(next);
        statezz = M3(next);
        state12 = M4(next);
        state13 = M5(next);
        state21 = M6(next);
        state23 = M7(next);
        state31 = M8(next);
        state32 = M9(next);
        states = [statexx stateyy statezz];
        states_shear = [state12 state13 state21 state23 state31 state32];
        color = rand(1,3);
        for i = 1:peaks
            dpx = statexx(i) - startin_ptx;
            dpy = stateyy(i) - startin_pty;
            dpz = statezz(i) - startin_ptz;
            set(0, 'CurrentFigure', f1)
            quiver3(startin_ptx,startin_pty,startin_ptz,dpx,dpy,dpz,0,'LineWidth',0.8,'Color', color);
            startin_ptx = startin_ptx + dpx;
            startin_pty = startin_pty + dpy;
            startin_ptz = startin_ptz + dpz;
        end
        
        trueStrain = {};
        trueStrain_shear = {};
        b = 1;

        %Here is the deformation vector for F11, F22 and F33 created:
        while b == 1
        et = 1e-5 + (1e-3-1e-5).*rand(1);
        dt = round(0.05 + (5-0.05).*rand(1),3);
        de1 = 1e-6 + (1e-4-1e-6).*rand(1);
        de2 = 1e-6 + (1e-4-1e-6).*rand(1); 
        for i = 1:3
            k = 1;
            for j = 2:peaks+1
                prev_peak = lpd(i,j-1);
                    current_peak = states(j-1,i);
                    de = de1;
                    if current_peak > prev_peak
                        current_strain = prev_peak:de:current_peak;
                    elseif current_peak < prev_peak
                        current_strain = prev_peak:-de:current_peak;
                    end
                trueStrain{i,k} = current_strain; 
                k = k +1;
                lpd(i,j) = current_peak;
            end
        end
        de1
        %Here is the deformation vector for F12,F13,F21... created:
        for i = 1:6
            k = 1;
            for j = 2:peaks+1
                prev_peak = lps(i,j-1);
                    current_peak = states_shear(j-1,i);
                    de = de2;
                    if current_peak > prev_peak
                        current_strain = prev_peak:de:current_peak;
                    elseif current_peak < prev_peak
                        current_strain = prev_peak:-de:current_peak;
                    end
                trueStrain_shear{i,k} = current_strain; 
                k = k +1;
                lps(i,j) = current_peak;
            end
        end
        de2
        %Collect strain vectors for F11,F22 and F33
        trueStrain1 ={};
        trueStrain2 = {};
        trueStrain3 = {};
        for i = 1:3
            for j = 1:size(trueStrain,2)
                switch i
                    case 1
                        trueStrain1{1,j} = trueStrain{i,j};
                    case 2
                        trueStrain2{1,j} = trueStrain{i,j};
                    case 3
                        trueStrain3{1,j} = trueStrain{i,j};
                end
            end
        end
        trueStrain11 = cell2mat(trueStrain1);
        trueStrain22 = cell2mat(trueStrain2);
        trueStrain33 = cell2mat(trueStrain3);
        %Collect strain vectors for F12,F13,F21 etc.
        trueStrain12 ={};
        trueStrain13 = {};
        trueStrain21 = {};
        trueStrain23 = {};
        trueStrain31 = {};
        trueStrain32 = {};
        for i = 1:6
            for j = 1:size(trueStrain_shear,2)
                switch i
                    case 1
                        trueStrain12{1,j} = trueStrain_shear{i,j};
                        trueStrain{4,j} = trueStrain_shear{i,j};
                    case 2
                        trueStrain13{1,j} = trueStrain_shear{i,j};
                        trueStrain{5,j} = trueStrain_shear{i,j};
                    case 3
                        trueStrain21{1,j} = trueStrain_shear{i,j};
                        trueStrain{6,j} = trueStrain_shear{i,j};
                    case 4
                        trueStrain23{1,j} = trueStrain_shear{i,j};
                        trueStrain{7,j} = trueStrain_shear{i,j};
                    case 5
                        trueStrain31{1,j} = trueStrain_shear{i,j};
                        trueStrain{8,j} = trueStrain_shear{i,j};
                    case 6
                        trueStrain32{1,j} = trueStrain_shear{i,j};
                        trueStrain{9,j} = trueStrain_shear{i,j};
                end
            end
        end
        trueStrain12 =cell2mat(trueStrain12);
        trueStrain13 = cell2mat(trueStrain13);
        trueStrain21 = cell2mat(trueStrain21);
        trueStrain23 = cell2mat(trueStrain23);
        trueStrain31 = cell2mat(trueStrain31);
        trueStrain32 = cell2mat(trueStrain32);

        % Create deformation gradient here just for one step to calculate
        % the strain and strain rate and check if we are within a defined
        % strain rate !!
        % If not continue the while loop again until we are within the
        % defined strain rate
        switch combination
            case 1
               F_diag = [1+de1 0 0;0 1 0;0 0 1]; 
            case 2
                F_diag = [1 0 0;0 1+de1 0;0 0 1];
            case 3
                F_diag = [1 0 0;0 1 0;0 0 1+de1];
            case 12 
                F_diag = [1+de1 0 0;0 1+de1 0;0 0 1];
            case 13
                F_diag = [1+de1 0 0;0 1 0;0 0 1+de1];
            case 23
                F_diag = [1 0 0;0 1+de1 0;0 0 1+de1];
            case 123
                F_diag = [1+de1 0 0;0 1+de1 0;0 0 1+de1];
            case 9
                F_diag = [1+de1 de2 de2;de2 1+de1 de2;de2 de2 1+de1];
            case 112
                F_diag = [1+de1 de2 0;0 1 0;0 0 1];
            case 113
                F_diag = [1+de1 0 de2;0 1 0;0 0 1];
            case 121
                F_diag = [1+de1 0 0;de2 1 0;0 0 1];
            case 1230
                F_diag = [1+de1 0 0;0 1 de2;0 0 1];
            case 131
                F_diag = [1+de1 0 0;0 1 0;de2 0 1];
            case 132
                F_diag = [1+de1 0 0;0 1 0;0 de2 1];
            case 212
                F_diag = [1 de2 0;0 1+de2 0;0 0 1];
            case 213
                F_diag = [1 0 de2;0 1+de1 0;0 0 1];
            case 221
                F_diag = [1 0 0;de2 1+de1 0;0 0 1];
            case 223
                F_diag = [1 0 0;0 1+de1 de2;0 0 1];
            case 231
                F_diag = [1 0 0;0 1+de1 0;de2 0 1];
            case 232
                F_diag = [1 0 0;0 1+de1 0;0 de2 1];
            case 312
                F_diag = [1 de2 0;0 1 0;0 0 1+de1];
            case 313
                F_diag = [1 0 de2;0 1 0;0 0 1+de1];
            case 321
                F_diag = [1 0 0;de2 1 0;0 0 1+de1];
            case 323
                F_diag = [1 0 0;0 1 de2;0 0 1+de1];
            case 331
                F_diag = [1 0 0;0 1 0;de2 0 1+de1];
            case 332
                F_diag = [1 0 0;0 1 0;0 de2 1+de1];                
        end
        % Calculate strain
        strain_diag = 1/2*(eye(3) - inv(F_diag*F_diag'));
        norm_s = norm(strain_diag,'fro');
        %Check if within strain rate
        if abs(norm_s/dt) > 1e-5 & abs(norm_s/dt) < 1e-3
            b = 0;
        else
            b = 1;
        end
        end

        % Make sure all deformation vectors have the same length

        ind = randi([1,9],1);
        trueStrain_copy = trueStrain;
        trueStrain_copy2 = trueStrain;
        for sim_types = 1
            if sim_types == 1
                switch ind;
                    case {1,4,5}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} =  cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {2,6,7}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} = cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {3,8,9}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} = cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                end
                trueStrain11 = trueStrain_final{1};
                trueStrain22 = trueStrain_final{2};
                trueStrain33 = trueStrain_final{3};
                trueStrain12 = trueStrain_final{4};
                trueStrain13 = trueStrain_final{5};
                trueStrain21 = trueStrain_final{6};
                trueStrain23 = trueStrain_final{7};
                trueStrain31 = trueStrain_final{8};
                trueStrain32 = trueStrain_final{9};
            else
                switch ind;
                    case {1,4,5}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {2,6,7}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {3,8,9}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                end
                trueStrain11 = trueStrain_final2{1};
                trueStrain22 = trueStrain_final2{2};
                trueStrain33 = trueStrain_final2{3};
                trueStrain12 = trueStrain_final2{4};
                trueStrain13 = trueStrain_final2{5};
                trueStrain21 = trueStrain_final2{6};
                trueStrain23 = trueStrain_final2{7};
                trueStrain31 = trueStrain_final2{8};
                trueStrain32 = trueStrain_final2{9};
            end
            %Run simulation
    %         for tempsat = 1:2

                for zitamat = 1
                    for wnpmat = 1
                        zita = zitam(zitamat);
                        wnp   = wnp_m(wnpmat);
    %                     Temper = temps(tempsat)
                        alphaZ = 1 + 0.057*zita.^2-9.5.*zita;
                        alphaT = 2 - exp(0.0126*(Temper - Tref));
                        params_F = [Temper Tref alphaT zita alphaZ];
                        wgf   = 0.0;              % glass fiber weight fraction
                        wp    = 1 - wnp - wgf;    % polymer weight fraction
                        % densities
                        ro_p  = 1.20;               % density of polymer (g/ml)
                        ro_np = 3.00;               % density of nanoparticle (g/ml)
                        ro_gf = 2.55;               % density of glass fiber (g/ml)
                        % nanoparticle volume fraction
                        vnp    = wnp * ro_p / (ro_np + wnp*ro_p - ro_np*wnp);                
                        if combination == 9
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt_mono,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    trueStrain12,trueStrain13,trueStrain21,trueStrain23,trueStrain31,trueStrain32,0);
                            toc  
                        elseif combination == 1
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc    
                        elseif combination == 2
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                        elseif combination == 3
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                        elseif combination == 12
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 13
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 23
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt_mono,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                        elseif combination == 123
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc
                        elseif combination == 112
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 113
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                        elseif combination == 121
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc             
                        elseif combination == 1230
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc         
                        elseif combination == 131
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                        elseif combination == 132
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc  
                
                        elseif combination == 212
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 213
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                        elseif combination == 221
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc             
                        elseif combination == 223
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc         
                        elseif combination == 231
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                        elseif combination == 232
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc
                
                        elseif combination == 312
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 313
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                        elseif combination == 321
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc             
                        elseif combination == 323
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc         
                        elseif combination == 331
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                        elseif combination == 332
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc            
                        end
                        sig_rheo = cell2mat(sigma);
                        F_rheo = cell2mat(F);
                        F_v_rheo = cell2mat(F_v);
                        F_vp_rheo = cell2mat(F_vp);
                        strain_rheo = cell2mat(strain);
                        b_mat = cell2mat(b_almansi);
                        time = t_tot';
                        C_rheo = cell2mat(C_pert)';
                        if length(sig_rheo) ~= length(strain_rheo)
                            continue
                        end 
                        if length(sig_rheo) < 1
                            continue
                        end                    
                        % Plot stress-strain 
                        set(0, 'CurrentFigure', f2)
                        hold on;
                        plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                        'LineWidth',1,...
                        'LineStyle','-','Color','b'); 
                        plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                        'LineWidth',1.5,...
                        'LineStyle','--','Color','r');
                        plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                        'LineWidth',2,...
                        'LineStyle','-.','Color','k');
                        
                
                        set(0,'CurrentFigure',f4)
                        hold on;
                        plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                        'LineWidth',1,...
                        'LineStyle','-','Color','b'); 
                        plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                        'LineWidth',1.5,...
                        'LineStyle','--','Color','r');
                        plot(strain_rheo(2:3:end,1),sig_rheo(2:3:end,1),'DisplayName','S21',...
                        'LineWidth',2,...
                        'LineStyle','-.','Color','k');    
                
                        plot(strain_rheo(2:3:end,3),sig_rheo(2:3:end,3),'DisplayName','S23',...
                        'LineWidth',1,...
                        'LineStyle','-','Color','g'); 
                        plot(strain_rheo(3:3:end,1),sig_rheo(3:3:end,1),'DisplayName','S31',...
                        'LineWidth',1.5,...
                        'LineStyle','--','Color','r');
                        plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S32',...
                        'LineWidth',2,...
                        'LineStyle','-.','Color','k'); 
                        
                        set(0, 'CurrentFigure', f3)
                        hold on;
                        for i = 1:3
                        switch i
                            case 1
                                plot(trueStrain11,'Color','b','LineStyle','-','LineWidth',1);
                            case 2
                                plot(trueStrain22,'Color','r','LineStyle','--','LineWidth',1.5);
                            case 3
                                plot(trueStrain33,'Color','k','LineStyle','-.','LineWidth',2);
                        end
                        end
                    %     pause(2);
                        
                        %Create input and output for stress-b mapping
                        inputs = {};
                        outputs = {};
                        k = 1;
                        for i = 1:length(F)
                                F1 = F{i,1};
                                b_1 = b_almansi{i,1};
                                current_F_1 = wrapper(F{i,1})';
                                current_stress_1 = wrapper_symm(sigma{i,1})';
                                current_b = wrapper_symm(b_1)';
                                inputs{1,k} = [current_b;dt;wnp;zita;Temper];
                                outputs{1,k} = [current_stress_1];
                                k = k + 1;
                        end
                        x = cell2mat(inputs);
                        t = cell2mat(outputs);
                        if any(t > abs(400), 'all')
                            sims
                            continue
                        else
                            target = sprintf('simulations_const_random_%d/input_target_ml_et_%d_%d_%d_%d.mat',peaks,sims,wnpmat,zitamat,sim_types)
                            save(target,'x','t');
                        end
                    end
                end
        end
%         end
        %Create input and output for F and C_perturbation mapping
        inputs = {};
        outputs = {};
        k = 1;

        sim_strain{1,sims} = trueStrain11;
        sim_strain{2,sims} = trueStrain22;
        sim_strain{3,sims} = trueStrain33;
    end
    target = sprintf('/simulations_const_random_%d/space.fig',peaks);
    saveas(figure(f1),[pwd target]);
    target = sprintf('/simulations_const_random_%d/load_path.fig',peaks);
    saveas(figure(f3),[pwd target]); 
    target = sprintf('/simulations_const_random_%d/results.fig',peaks);
    saveas(figure(f2),[pwd target]);   
    %Validate the predicted stress strain from ML 
elseif validate == 1
    close all;
    load('lstm_3_seqtoseq_2150_64.mat');
    printout(net,muX,sigmaX,muT,sigmaT,64);
    try
        rmdir validation
    end
    mkdir validation
    for sims = 1:100
        close all;
%         combination =  randsample(pos_comb,1);
        combination =  randsample(pos_comb,1); %Randomly select combination of F  
        combination = 9;
%         combination = 1;
        %Starting points
        startin_ptx = 1.0; 
        startin_pty = 1.0;
        startin_ptz = 1.0;
        % Select randomly states for the simulation from the uniform points
        next = randi([1 length(X0)],peaks,1);
        statexx = M1(next);
        stateyy = M2(next);
        statezz = M3(next);
        state12 = M4(next);
        state13 = M5(next);
        state21 = M6(next);
        state23 = M7(next);
        state31 = M8(next);
        state32 = M9(next);
        states = [statexx stateyy statezz];
        states_shear = [state12 state13 state21 state23 state31 state32];
%         color = rand(1,3);
%         for i = 1:peaks
%             dpx = statexx(i) - startin_ptx;
%             dpy = stateyy(i) - startin_pty;
%             dpz = statezz(i) - startin_ptz;
%             quiver3(startin_ptx,startin_pty,startin_ptz,dpx,dpy,dpz,0,'LineWidth',0.8,'Color', color);
%             startin_ptx = startin_ptx + dpx;
%             startin_pty = startin_pty + dpy;
%             startin_ptz = startin_ptz + dpz;
%         end
        
        trueStrain = {};
        trueStrain_shear = {};
        b = 1;

        %Here is the deformation vector for F11, F22 and F33 created:
        %Here is the deformation vector for F11, F22 and F33 created:
        while b == 1
        et = 1e-5 + (1e-3-1e-5).*rand(1);
        dt = round(0.5 + (5-0.5).*rand(1),3);
        de1 = 8e-5 + (1e-4-8e-5).*rand(1);
        de2 = 5e-5 + (1e-4-5e-5).*rand(1); 
        for i = 1:3
            k = 1;
            for j = 2:peaks+1
                prev_peak = lpd(i,j-1);
                    current_peak = states(j-1,i);
                    de = de1;
                    if current_peak > prev_peak
                        current_strain = prev_peak:de:current_peak;
                    elseif current_peak < prev_peak
                        current_strain = prev_peak:-de:current_peak;
                    end
                trueStrain{i,k} = current_strain; 
                k = k +1;
                lpd(i,j) = current_peak;
            end
        end
        de1
        %Here is the deformation vector for F12,F13,F21... created:
        for i = 1:6
            k = 1;
            for j = 2:peaks+1
                prev_peak = lps(i,j-1);
                    current_peak = states_shear(j-1,i);
                    de = de2;
                    if current_peak > prev_peak
                        current_strain = prev_peak:de:current_peak;
                    elseif current_peak < prev_peak
                        current_strain = prev_peak:-de:current_peak;
                    end
                trueStrain_shear{i,k} = current_strain; 
                k = k +1;
                lps(i,j) = current_peak;
            end
        end
        de2
        %Collect strain vectors for F11,F22 and F33
        trueStrain1 ={};
        trueStrain2 = {};
        trueStrain3 = {};
        for i = 1:3
            for j = 1:size(trueStrain,2)
                switch i
                    case 1
                        trueStrain1{1,j} = trueStrain{i,j};
                    case 2
                        trueStrain2{1,j} = trueStrain{i,j};
                    case 3
                        trueStrain3{1,j} = trueStrain{i,j};
                end
            end
        end
        trueStrain11 = cell2mat(trueStrain1);
        trueStrain22 = cell2mat(trueStrain2);
        trueStrain33 = cell2mat(trueStrain3);
        %Collect strain vectors for F12,F13,F21 etc.
        trueStrain12 ={};
        trueStrain13 = {};
        trueStrain21 = {};
        trueStrain23 = {};
        trueStrain31 = {};
        trueStrain32 = {};
        for i = 1:6
            for j = 1:size(trueStrain_shear,2)
                switch i
                    case 1
                        trueStrain12{1,j} = trueStrain_shear{i,j};
                        trueStrain{4,j} = trueStrain_shear{i,j};
                    case 2
                        trueStrain13{1,j} = trueStrain_shear{i,j};
                        trueStrain{5,j} = trueStrain_shear{i,j};
                    case 3
                        trueStrain21{1,j} = trueStrain_shear{i,j};
                        trueStrain{6,j} = trueStrain_shear{i,j};
                    case 4
                        trueStrain23{1,j} = trueStrain_shear{i,j};
                        trueStrain{7,j} = trueStrain_shear{i,j};
                    case 5
                        trueStrain31{1,j} = trueStrain_shear{i,j};
                        trueStrain{8,j} = trueStrain_shear{i,j};
                    case 6
                        trueStrain32{1,j} = trueStrain_shear{i,j};
                        trueStrain{9,j} = trueStrain_shear{i,j};
                end
            end
        end
        trueStrain12 =cell2mat(trueStrain12);
        trueStrain13 = cell2mat(trueStrain13);
        trueStrain21 = cell2mat(trueStrain21);
        trueStrain23 = cell2mat(trueStrain23);
        trueStrain31 = cell2mat(trueStrain31);
        trueStrain32 = cell2mat(trueStrain32);

        % Create deformation gradient here just for one step to calculate
        % the strain and strain rate and check if we are within a defined
        % strain rate !!
        % If not continue the while loop again until we are within the
        % defined strain rate
        switch combination
            case 1
               F_diag = [1+de1 0 0;0 1 0;0 0 1]; 
            case 2
                F_diag = [1 0 0;0 1+de1 0;0 0 1];
            case 3
                F_diag = [1 0 0;0 1 0;0 0 1+de1];
            case 12 
                F_diag = [1+de1 0 0;0 1+de1 0;0 0 1];
            case 13
                F_diag = [1+de1 0 0;0 1 0;0 0 1+de1];
            case 23
                F_diag = [1 0 0;0 1+de1 0;0 0 1+de1];
            case 123
                F_diag = [1+de1 0 0;0 1+de1 0;0 0 1+de1];
            case 9
                F_diag = [1+de1 de2 de2;de2 1+de1 de2;de2 de2 1+de1];
            case 112
                F_diag = [1+de1 de2 0;0 1 0;0 0 1];
            case 113
                F_diag = [1+de1 0 de2;0 1 0;0 0 1];
            case 121
                F_diag = [1+de1 0 0;de2 1 0;0 0 1];
            case 1230
                F_diag = [1+de1 0 0;0 1 de2;0 0 1];
            case 131
                F_diag = [1+de1 0 0;0 1 0;de2 0 1];
            case 132
                F_diag = [1+de1 0 0;0 1 0;0 de2 1];
            case 212
                F_diag = [1 de2 0;0 1+de2 0;0 0 1];
            case 213
                F_diag = [1 0 de2;0 1+de1 0;0 0 1];
            case 221
                F_diag = [1 0 0;de2 1+de1 0;0 0 1];
            case 223
                F_diag = [1 0 0;0 1+de1 de2;0 0 1];
            case 231
                F_diag = [1 0 0;0 1+de1 0;de2 0 1];
            case 232
                F_diag = [1 0 0;0 1+de1 0;0 de2 1];
            case 312
                F_diag = [1 de2 0;0 1 0;0 0 1+de1];
            case 313
                F_diag = [1 0 de2;0 1 0;0 0 1+de1];
            case 321
                F_diag = [1 0 0;de2 1 0;0 0 1+de1];
            case 323
                F_diag = [1 0 0;0 1 de2;0 0 1+de1];
            case 331
                F_diag = [1 0 0;0 1 0;de2 0 1+de1];
            case 332
                F_diag = [1 0 0;0 1 0;0 de2 1+de1];                
        end
        % Calculate strain
        strain_diag = 1/2*(eye(3) - inv(F_diag*F_diag'));
        norm_s = norm(strain_diag,'fro');
        %Check if within strain rate
        if abs(norm_s/dt) > 1e-5 & abs(norm_s/dt) < 1e-3
            b = 0;
        else
            b = 1;
        end
        end

        % Make sure all deformation vectors have the same length
        ind = randi([1,9],1);
        trueStrain_copy = trueStrain;
        trueStrain_copy2 = trueStrain;
        for sim_types = 1:2
            if sim_types == 1
                switch ind;
                    case {1,4,5}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} =  cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {2,6,7}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} = cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {3,8,9}
                        slength = cat(2,trueStrain_copy{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy{i,:}))
                                if i <= 3
                                    add_to = ones(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                else
                                    add_to = zeros(length(slength) - length(cat(2,trueStrain_copy{i,:})),1);
                                end
                                trueStrain_copy{i,1} = [add_to' trueStrain_copy{i,1}];
                                trueStrain_final{i} = cat(2,trueStrain_copy{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy{i,:});
                                trueStrain_final{i} = [current_s(1:end-diff_l)];
                            end
                        end
                end
                trueStrain11 = trueStrain_final{1};
                trueStrain22 = trueStrain_final{2};
                trueStrain33 = trueStrain_final{3};
                trueStrain12 = trueStrain_final{4};
                trueStrain13 = trueStrain_final{5};
                trueStrain21 = trueStrain_final{6};
                trueStrain23 = trueStrain_final{7};
                trueStrain31 = trueStrain_final{8};
                trueStrain32 = trueStrain_final{9};
            else
                switch ind;
                    case {1,4,5}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {2,6,7}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                    case {3,8,9}
                        slength = cat(2,trueStrain_copy2{ind,:});
                        for i = 1:9
                            if length(slength) > length(cat(2,trueStrain_copy2{i,:}))
                                add_to = zeros(length(slength) - length(cat(2,trueStrain_copy2{i,:})),1);
                                add = cat(2,trueStrain_copy2{i,:});
                                add_to(:) = add(end);
                                trueStrain_copy2{i,end} = [trueStrain_copy2{i,end} add_to'];
                                trueStrain_final2{i} = cat(2,trueStrain_copy2{i,:});
                            else
                                diff_l = length(cat(2,trueStrain_copy2{i,:})) - length(slength);
                                current_s = cat(2,trueStrain_copy2{i,:});
                                trueStrain_final2{i} = [current_s(1:end-diff_l)];
                            end
                        end
                end
                trueStrain11 = trueStrain_final2{1};
                trueStrain22 = trueStrain_final2{2};
                trueStrain33 = trueStrain_final2{3};
                trueStrain12 = trueStrain_final2{4};
                trueStrain13 = trueStrain_final2{5};
                trueStrain21 = trueStrain_final2{6};
                trueStrain23 = trueStrain_final2{7};
                trueStrain31 = trueStrain_final2{8};
                trueStrain32 = trueStrain_final2{9};
            end
            %Run simulation
    %         for tempsat = 1:2
                sims_d = 1;
                for zitamat = 1:2
                    for wnpmat = 1:2
                        zita = zitam(zitamat);
                        wnp   = wnp_m(wnpmat);
    %                     Temper = temps(tempsat)
                        alphaZ = 1 + 0.057*zita.^2-9.5.*zita;
                        alphaT = 2 - exp(0.0126*(Temper - Tref));
                        params_F = [Temper Tref alphaT zita alphaZ wnp];
                        wgf   = 0.0;              % glass fiber weight fraction
                        wp    = 1 - wnp - wgf;    % polymer weight fraction
                        % densities
                        ro_p  = 1.20;               % density of polymer (g/ml)
                        ro_np = 3.00;               % density of nanoparticle (g/ml)
                        ro_gf = 2.55;               % density of glass fiber (g/ml)
                        % nanoparticle volume fraction
                        vnp    = wnp * ro_p / (ro_np + wnp*ro_p - ro_np*wnp);
                        wnp
                        zita
                        %Run simulation
                        if combination == 9
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    trueStrain12,trueStrain13,trueStrain21,trueStrain23,trueStrain31,trueStrain32,0);
                            toc
                            time_cons(sims,sims_d) = toc;
                            
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de_mono,d0,dt_mono,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    trueStrain12,trueStrain13,trueStrain21,trueStrain23,trueStrain31,trueStrain32,...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims,sims_d) = toc;
                            sims_d = sims_d + 1;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;
                        elseif combination == 0
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt_mono,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                        elseif combination == 1
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc    
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;              
                        elseif combination == 2
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;               
                        elseif combination == 3
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc   
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;               
                        elseif combination == 123
                            sprintf('Start of const. Model')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                             % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;
                        elseif combination == 13
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;             
                        elseif combination == 23
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;             
                        elseif combination == 12
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de_mono,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc     
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;
                        elseif combination == 112
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 113
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 121
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 1230
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc
                            time_cons(sims) = toc;
                            sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 131
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 132
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc  
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,trueStrain11,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;
                        elseif combination == 212
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 213
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 221
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 223
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 231
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 232
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),trueStrain22,ones(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;
                        elseif combination == 312
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    trueStrain12,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 313
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),trueStrain13,zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 321
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc    
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain21,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 323
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),0);
                            toc  
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    trueStrain23,zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 331
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),trueStrain31,zeros(1,length(trueStrain11)),...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        elseif combination == 332
                             tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,0);
                            toc 
                            time_cons(sims) = toc;
                           sig_rheo = cell2mat(sigma);
                            F_rheo = cell2mat(F);
                            F_v_rheo = cell2mat(F_v);
                            F_vp_rheo = cell2mat(F_vp);
                            strain_rheo = cell2mat(strain);
                            b_mat = cell2mat(b_almansi);
                            time = t_tot';   
                            sprintf('Start of ML')
                            tic
                            [stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
                                    strain,trueStrain,b_almansi] = model_vevpd_fmix_ml(params,de,d0,dt,vnp,params_F,...
                                    nnetwork,0,et_load,et_unload,ones(1,length(trueStrain11)),ones(1,length(trueStrain11)),trueStrain33,...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),...
                                    zeros(1,length(trueStrain11)),zeros(1,length(trueStrain11)),trueStrain32,...
                                    muX,sigmaX,muT,sigmaT,net,rdn);              
                            toc
                            time_ml(sims) = toc;
                            sig_ml = cell2mat(sigma);
                            F_ml = cell2mat(F);
                            strain_ml = cell2mat(strain);
                            % Plot stress-strain 
                            h = figure('Name','Rheology model vs Machine Learning');
                            tiledlayout(2,6)
                            ax1 = nexttile;
                            title(ax1,'S11')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,1),sig_rheo(1:3:end,1),'DisplayName','S11',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b'); 
                            ax1 = nexttile;
                            title(ax1,'S22')
                            ax1.FontSize = 14;
                            plot(strain_rheo(2:3:end,2),sig_rheo(2:3:end,2),'DisplayName','S22',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');      
                            ax1 = nexttile;
                            title(ax1,'S33')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,3),sig_rheo(3:3:end,3),'DisplayName','S33',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');   
                            ax1 = nexttile;
                            title(ax1,'S12')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,2),sig_rheo(1:3:end,2),'DisplayName','S12',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S13')
                            ax1.FontSize = 14;
                            plot(strain_rheo(1:3:end,3),sig_rheo(1:3:end,3),'DisplayName','S13',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                            ax1 = nexttile;
                            title(ax1,'S23')
                            ax1.FontSize = 14;
                            plot(strain_rheo(3:3:end,2),sig_rheo(3:3:end,2),'DisplayName','S23',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','b');
                
                            ax1 = nexttile;
                            title(ax1,'S11 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,1),sig_ml(1:3:end,1),'DisplayName','S11 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r'); 
                            ax1 = nexttile;
                            title(ax1,'S22 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(2:3:end,2),sig_ml(2:3:end,2),'DisplayName','S22 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');      
                            ax1 = nexttile;
                            title(ax1,'S33 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,3),sig_ml(3:3:end,3),'DisplayName','S33 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');  
                            ax1 = nexttile;
                            title(ax1,'S12 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,2),sig_ml(1:3:end,2),'DisplayName','S12 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S13 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(1:3:end,3),sig_ml(1:3:end,3),'DisplayName','S13 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                            ax1 = nexttile;
                            title(ax1,'S23 - ML')
                            ax1.FontSize = 14;
                            plot(strain_ml(3:3:end,2),sig_ml(3:3:end,2),'DisplayName','S23 - ML',...
                            'LineWidth',1,...
                            'LineStyle','-','Color','r');
                
                            target = sprintf('/validation/sim_%d.fig',sims);
                            saveas(figure(h),[pwd target]);
                            hold off;            
                        end
%                         close all;
                        rheo_sv = reshape(sig_rheo,1,[]);
                        ml_sv = reshape(sig_ml,1,[]);
                        if length(rheo_sv) ~= length(ml_sv)
                            continue
                        end
                        Err = rheo_sv - ml_sv;      
                        me_c{sims} =  mean(abs(Err));
                        RMSE_c{sims}= sqrt(mean(Err.^2));
                        std_c{sims} = std(Err);
                    end
                end
        end

        
%         set(0, 'CurrentFigure', f3)
%         hold on;
%         for i = 1:3
%         switch i
%             case 1
%                 plot(trueStrain11,'Color','b');
%             case 2
%                 plot(trueStrain22,'Color','r');
%             case 3
%                 plot(trueStrain33);
%         end
%         end  
%         pause(2)
%         close(h);
    end
    index = cellfun(@isempty, RMSE_c) == 0;
    RMSE_c = RMSE_c(index);
    me_c = me_c(index);
    rmsev = cell2mat(RMSE_c)';
    mev = cell2mat(me_c)';
    close all;
    h = figure('Name','Rheology model vs Machine Learning - Dist_RMSE');
    histogram(rmsev)
    hold on;
%     ylim([0 1.1]);
    target = sprintf('/validation/dist_rmse.fig');
    saveas(figure(h),[pwd target]);

    h = figure('Name','Rheology model vs Machine Learning - Dist_MAE');
    histogram(mev,10)
    hold on;
    grid on;
    xlabel('Data Value', 'FontSize', 15);
    ylabel('Count', 'FontSize', 15);
    % Compute mean and standard deviation.
    mu = mean(mev)
    sigma = std(mev)
    % Indicate those on the plot.
    xline(mu, 'Color', 'g', 'LineWidth', 2);
    xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    ylim([0, sims]); % Give some headroom above the bars.
    yl = ylim;
    sMean = sprintf('  Mean = %.3f\n  SD = %.3f', mu, sigma);
    % Position the text 90% of the way from bottom to top.
    text(mu, 0.9*yl(2), sMean, 'Color', 'r', ...
	    'FontWeight', 'bold', 'FontSize', 12, ...
	    'EdgeColor', 'b');
    sMean2= sprintf('Histogram for MAE.  Mean = %.3f.  SD = %.3f', mu, sigma);
    title(sMean2, 'FontSize', 15);

%     ylim([0 1.1]);
    target = sprintf('/validation/dist_mae.fig');
    saveas(figure(h),[pwd target]);
    fullfile_tex = ['/validation/dist_mae' '.tex'];
    matlab2tikz(fullfile_tex);
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

function [y] = wrapper_c(matrix)
    C1111 = matrix(1,1);
    C2222 = matrix(5,5);
    C3333 = matrix(9,9);
    C1212 = matrix(4,1);
    C1313 = matrix(3,7);
    C2323 = matrix(6,8);
    C1122 = matrix(2,2);
    C1133 = matrix(3,3);
    C2233 = matrix(6,6);
    y = [C1111;C1122;C1133;C3333;C2323;C1212;C1313;C2222;C2233];
end

function [y] = wrapper_symm(matrix)
    S11 = matrix(1,1);
    S12 = matrix(1,2);
    S13 = matrix(1,3);
    S22 = matrix(2,2);
    S23 = matrix(2,3);
    S33 = matrix(3,3);

    y = [S11 S12 S13 S22 S23 S33];
end
