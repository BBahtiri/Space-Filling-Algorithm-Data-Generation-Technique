function [stress,F_tot,F_ve_t,F_v_t,F_e_t,F_vp_t,J_tot,timeVec,...
    T_model_all,T_model_networkA,T_model_networkB,T_model_volumetric,d_total,...
    strain_gl,trueStrain1,b_alm] = model_calibr(params,de,d0,dt,vp,params_F,nnetwork,...
    cyclic_loading,et_load,et_unload,trueStrain1,trueStrain2,trueStrain3,trueStrain12,...
    trueStrain13,trueStrain21,trueStrain23,trueStrain31,trueStrain32,...
    muX,sigmaX,muT,sigmaT,net,rdn)
format longG
global eps0
eps0 = 0;
lc_max = 1;
d_t    = d0;
T      = zeros(3);
F22    = single(1); % Set F22 to 1 
T_model_loading = zeros(length(trueStrain1),1);
T_model_networkA = {};
T_model_networkB = {};
T_model_volumetric = {};
T_model_all = {};
zita = params_F(4);
wnp = params_F(6);
strain_gl = {};
d_total         = zeros(length(trueStrain1),1);
J_tot         = zeros(length(trueStrain1),1);
t_tot = zeros(length(trueStrain1),1);
F_tot = {};
F_ve_t = {};
F_v_t = {};
F_e_t = {};
F_vp_t = {};
timeVec = zeros(length(trueStrain1),1);

TB_t = zeros(3);
strain = 0;

YPred_or = {};
% net = resetState(net);
% net0 = net;
% F_init = ([1+trueStrain(1) 1+trueStrain(2) 1+trueStrain(3) 1+trueStrain(4) 1+trueStrain(5) ...
%      1+trueStrain(6) 1+trueStrain(7) 1+trueStrain(8) 1+trueStrain(9)]...
%       -muX)./sigmaX;
% %     1+trueStrain(5) 1+trueStrain(6) 1+trueStrain(7) 1+trueStrain(8) 1+trueStrain(9)...
% %     1+trueStrain(10) 1+trueStrain(11) 1+trueStrain(12) 1+trueStrain(13)...
% %     1+trueStrain(14) 1+trueStrain(15) 1+trueStrain(16)];
% F_init = 1.0;
% [net,~] = predictAndUpdateState(net,F_init);

%% Get the weights and bias of the fully connected layer
% Get the weights and bias of the LSTM layers

if size(net.Layers,1) == 5
    %Get weights for first LSTM
    iweights_lstm1 = net.Layers(2, 1).InputWeights  ;
    rweights_lstm1 = net.Layers(2, 1).RecurrentWeights  ;
    bias_lstm1 = net.Layers(2, 1).Bias  ;
    initial_cstate1 = net.Layers(2, 1).CellState  ;
    initial_hstate1 = net.Layers(2, 1).HiddenState  ;
    % Initiate cell state and hidden state:
    cstate1 = initial_cstate1;
    hstate1 = initial_hstate1;
    %Get weights for second LSTM
    iweights_lstm2 = net.Layers(3, 1).InputWeights  ;
    rweights_lstm2 = net.Layers(3, 1).RecurrentWeights  ;
    bias_lstm2 = net.Layers(3, 1).Bias  ;
    initial_cstate2 = net.Layers(3, 1).CellState  ;
    initial_hstate2 = net.Layers(3, 1).HiddenState  ;
    % Initiate cell state and hidden state:
    cstate2 = initial_cstate2;
    hstate2 = initial_hstate2;
    weights_cl = net.Layers(4, 1).Weights;
    bias_cl = net.Layers(4, 1).Bias;
elseif size(net.Layers,1) == 4
    %Get weights for first LSTM
    iweights_lstm = net.Layers(2, 1).InputWeights  ;
    rweights_lstm = net.Layers(2, 1).RecurrentWeights  ;
    bias_lstm = net.Layers(2, 1).Bias  ;
    initial_cstate = net.Layers(2, 1).CellState  ;
    initial_hstate = net.Layers(2, 1).HiddenState  ;    
    % Initiate cell state and hidden state:
    cstate = initial_cstate;
    hstate = initial_hstate;
    weights_cl = net.Layers(3, 1).Weights;
    bias_cl = net.Layers(3, 1).Bias;
elseif size(net.Layers,1) == 6
    %Get weights for first LSTM
    iweights_lstm1 = net.Layers(2, 1).InputWeights  ;
    rweights_lstm1 = net.Layers(2, 1).RecurrentWeights  ;
    bias_lstm1 = net.Layers(2, 1).Bias  ;
    initial_cstate1 = net.Layers(2, 1).CellState  ;
    initial_hstate1 = net.Layers(2, 1).HiddenState  ;
    % Initiate cell state and hidden state:
    cstate1 = initial_cstate1;
    hstate1 = initial_hstate1;
    %Get weights for second LSTM
    iweights_lstm2 = net.Layers(3, 1).InputWeights  ;
    rweights_lstm2 = net.Layers(3, 1).RecurrentWeights  ;
    bias_lstm2 = net.Layers(3, 1).Bias  ;
    initial_cstate2 = net.Layers(3, 1).CellState  ;
    initial_hstate2 = net.Layers(3, 1).HiddenState  ;
    % Initiate cell state and hidden state:
    cstate2 = initial_cstate2;
    hstate2 = initial_hstate2;
    %Get weights for third  LSTM
    iweights_lstm3 = net.Layers(4, 1).InputWeights  ;
    rweights_lstm3 = net.Layers(4, 1).RecurrentWeights  ;
    bias_lstm3 = net.Layers(4, 1).Bias  ;
    initial_cstate3 = net.Layers(4, 1).CellState  ;
    initial_hstate3 = net.Layers(4, 1).HiddenState  ;   
    % Initial cell state and hidden state
    cstate3 = initial_cstate3;
    hstate3 = initial_hstate3;    
    weights_cl = net.Layers(5, 1).Weights;
    bias_cl = net.Layers(5, 1).Bias;    
end

% epsilon = 1e-4;
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

%%

muX = single(muX);
sigmaX = single(sigmaX);
muT = single(muT);
sigmaT = single(sigmaT);
k = 1;
if cyclic_loading == 0
for j = 1:length(trueStrain1)
%         F11 = 1 + trueStrain(j); % Increment for F11
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
        b = wrapper_symm(F*F');
        b_alm = F*F';
        wF = wrapper(F);
%         inputs_vec = [([b(1);b(4);b(6)]-muX)./sigmaX ];
        inputs_vec = [([b(1);b(2);b(3);b(4);b(5);b(6);dt;wnp;zita]-muX)./sigmaX ];
        inputs = single(inputs_vec);
        J  = det(F);
        Fb = J^(-1/3) * F;
        wF = wrapper(F);
        Cb = F'*F;
        E_gl = 1/2*(Cb-eye(3));
        
        if size(net.Layers,1) == 5
            [Y, C, G] = lstm_forward(inputs,iweights_lstm1,rweights_lstm1,...
                    bias_lstm1,cstate1,hstate1);
            cstate1 = C;
            hstate1 = Y;
            [Y, C, G] = lstm_forward(Y,iweights_lstm2,rweights_lstm2,...
                    bias_lstm2,cstate2,hstate2);
            cstate2 = C;
            hstate2 = Y;    
            YPred_own = weights_cl*Y+bias_cl;
        elseif size(net.Layers,1) == 4
            [Y, C, G] = lstm_forward(inputs,iweights_lstm,rweights_lstm,...
                    bias_lstm,cstate,hstate);
            cstate = C;
            hstate = Y;
            YPred_own = weights_cl*Y+bias_cl;
        elseif size(net.Layers,1) == 6
            [Y, C, G] = lstm_forward(inputs,iweights_lstm1,rweights_lstm1,...
                    bias_lstm1,cstate1,hstate1);
            cstate1 = C;
            hstate1 = Y;
            [Y, C, G] = lstm_forward(Y,iweights_lstm2,rweights_lstm2,...
                    bias_lstm2,cstate2,hstate2);
            cstate2 = C;
            hstate2 = Y;
            [Y, C, G] = lstm_forward(Y,iweights_lstm3,rweights_lstm3,...
            bias_lstm3,cstate3,hstate3);
            cstate3 = C;
            hstate3 = Y;    
            YPred_own = weights_cl*Y+bias_cl;            
        end
%         [net,YPred_own] = predictAndUpdateState(net,inputs);
        YPred_or{1,1} = YPred_own.*sigmaT + muT;
        Td = cell2mat(YPred_or);
        Td = symmetrize(Td);
        T = (1 - d_t) * Td;
%          if (lc_t > lc_max)
%              lc_max = lc_t;
%             d_t    = 1 - exp(Ad * (1 - lc_max) );        
%          end
        Tdev                 = Dev(T);
        T_model_loading(k) = T(1,1);

        T_model_all{k,1} = T;
        strain_gl{k,1} = E_gl;
        timeVec(k+1) = timeVec(j)+dt;
        F_tot{k,1 } = F;
        J_tot(k) = J;
        stress = T_model_loading;
        k = k +1;
end
elseif cyclic_loading == 1
    strain_cycles = [0.0087 0.012 0.0155 ...
    0.0192 0.022 0.027 0.032];
%     strain_cycles = [0.0025 0.005 0.0075 ...
%     0.01 0.0125 0.015 0.0175];
%     strain_cycles = [0.01 0.02 0.03 ...
%     0.04 0.05 0.06 0.07];
    trueStrain = [];
    cycle = 1;
    cycles = 7;
    loading = 1;
    unloading = 0;
    timeVec = [0, dt];
    mul = 5;
    de = single(mul*5e-5);
    de2 = single(mul*6e-5);
    de3 = single(mul*3e-5);
    dt = de/et_load;
    trueStrain = [0];
    trueStrain2 = [0];
    trueStrain3 = [0];
    stress = [0];
    run = 1;
    j = 1;
    
while run == 1
    if loading == 1
        timeVec = [0, dt];
        time0 = timeVec(end-1);
        time1 = timeVec(end);
        F11 = single(1.0) + trueStrain(end);
%         F22 = single(1.0) + trueStrain2(end);
%         F33 = single(1.0) + trueStrain3(end);
        F22 = 1.0;
        F33 = 1.0;
        
        ep  = 1e-4;
        iter_stress = 0; 

                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;
                
                F = [F11 0 0; 0 F22 0; 0 0 F33]; % Incremental F at t+1
                b_alm = F*F';
                b = wrapper_symm(F*F');
                wF = wrapper(F);
                inputs_vec = [([b(1)]-muX)./sigmaX ];
                inputs = inputs_vec;
                J  = det(F);
                Fb = J^(-1/3) * F;  
                Cb = F'*F;
                E_gl = 1/2*(Cb-eye(3));  
                
%                 [Y, C, G] = lstm_forward(inputs,iweights_lstm,rweights_lstm,...
%                             bias_lstm,cstate,hstate);
%                 cstate = C;
%                 hstate = Y;
%                 YPred_own = weights_cl*Y+bias_cl;
                
                  [net,YPred_own] = predictAndUpdateState(net,inputs);
                
                YPred_or{1,1} = YPred_own.*sigmaT + muT;
                Td = [YPred_or{1,1}(1,1);0;0; YPred_or{1,1}(2,1);0; YPred_or{1,1}(3,1)];
                Td = symmetrize(Td);
                
                
                T = (1 - d_t) * Td;


%          if (lc_t > lc_max)
%              lc_max = lc_t;
%              d_t    = 1 - exp(Ad * (1 - lc_max) );        
%          end
            stress = [stress T(1,1)];
            T_model_all{j,1} = Td;
 
            strain_gl{j,1} = E_gl;
            %timeVec(j+1) = timeVec(j)+dt;
            F_tot{j,1 } = F;

            J_tot(j) = J;
            j = j + 1;
         if trueStrain(end) > strain_cycles(cycle)
             loading = 0;
             unloading = 1;
             timeVec = [0, dt];
             de = single(mul*5e-5);
             de2 = single(mul*6e-5);
             de3 = single(mul*3e-5);
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
        F11 = single(1.0) + trueStrain(end);
%         F22 = single(1.0) + trueStrain2(end);
%         F33 = single(1.0) + trueStrain3(end);
        F22 = 1.0;
        F33 = 1.0;
        ep  = 1e-4;
        iter_stress = 0; 
                Tn = T(2,2); % Stress in uniaxial direction
                iter_stress = iter_stress + 1;

                F = [F11 0 0; 0 F22 0; 0 0 F33]; % Incremental F at t+1
                b_alm = F*F';
                b = wrapper_symm(F*F');
                wF = wrapper(F);
                inputs_vec = [([b(1)]-muX)./sigmaX ];
                inputs = inputs_vec;
                J  = det(F);
                Fb = J^(-1/3) * F;  
                Cb = F'*F;
                E_gl = 1/2*(Cb-eye(3));                
                
%                 [Y, C, G] = lstm_forward(inputs,iweights_lstm,rweights_lstm,...
%                             bias_lstm,cstate,hstate);
%                 cstate = C;
%                 hstate = Y;
%                 YPred_own = weights_cl*Y+bias_cl;
                
                 [net,YPred_own] = predictAndUpdateState(net,inputs);
                
                YPred_or{1,1} = YPred_own.*sigmaT + muT;
                Td = [YPred_or{1,1}(1,1);0;0; YPred_or{1,1}(2,1);0; YPred_or{1,1}(3,1)];
                Td = symmetrize(Td);
                
                
                T = (1 - d_t) * Td;                


%          if (lc_t > lc_max)
%              lc_max = lc_t;
%              d_t    = 1 - exp(Ad * (1 - lc_max) );        
%          end
            stress = [stress T(1,1)];
            T_model_all{j,1} = Td;

            strain_gl{j,1} = E_gl;
            %timeVec(j+1) = timeVec(j)+dt;
            F_tot{j,1 } = F;
            J_tot(j) = J;
            j = j + 1;
         if stress(end) < 1e-4
             loading = 1;
             unloading = 0;
             de = single(mul*5e-5);
             de2 = single(mul*6e-5);
             de3 = single(mul*3e-5);
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
    y = [F11 F12 F13 F21 F22 F23 F31 F32 F33]';
end

function [y] = wrapp_back(vector)
    F11 = vector(1);
    F12 = vector(2);
    F13 = vector(3);
    F21 = vector(4);
    F22 = vector(5);
    F23 = vector(6);
    F31 = vector(7);
    F32 = vector(8);
    F33 = vector(9);
    y = [F11 F12 F13; F21 F22 F23; F31 F32 F33];
end

function [y] = symmetrize(vector)
    s11 = vector(1);
    s12 = vector(2);
    s13 = vector(3);
    s22 = vector(4);
    s23 = vector(5);
    s33 = vector(6);
    y = [s11 s12 s13; s12 s22 s23; s13 s23 s33];
end

function [y] = normalize_input(vector,max_x,min_x)
    x_new = zeros(0,length(vector));
    for i = 1:length(vector)
        current_inp = vector(i);
        max_x_i = max_x(i);
        min_x_i = min_x(i);

        new_inp = 0.1 + 0.8*((current_inp - min_x_i)/(max_x_i - min_x_i));
        
        x_new(i) = new_inp;

    end
    y = x_new';
end

function [y] = normalize_back(vector,max_t,min_t)
    x_new = zeros(0,length(vector));
    for i = 1:length(vector)
        new_out = vector(i);
        max_t_i = max_t(i);
        min_t_i = min_t(i);

%         new_out = 0.1 + 0.8*((current_out - min_t_i)/(max_t_i - min_t_i));
        
        current_out = (new_out/0.8)*(max_t_i - min_t_i) - 0.1 + min_t_i;
        
        x_new(i) = current_out;

    end
    y = x_new';
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

end

