function printout(net,muX,sigmaX,muT,sigmaT,n)
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
    networks = sprintf('network_2%d_%d',length(cstate1),n);
    mkdir(networks);
    cd (networks)
    writematrix(single(muX),'muX.txt','Delimiter',' ');
    writematrix(single(muT),'muT.txt','Delimiter',' ');
    writematrix(single(sigmaX),'sigmaX.txt','Delimiter',' ');
    writematrix(single(sigmaT),'sigmaT.txt','Delimiter',' ');
    writematrix(iweights_lstm1,'iweights_lstm1.txt','Delimiter',' ');
    writematrix(iweights_lstm2,'iweights_lstm2.txt','Delimiter',' ');
    writematrix(rweights_lstm1,'rweights_lstm1.txt','Delimiter',' ');
    writematrix(rweights_lstm2,'rweights_lstm2.txt','Delimiter',' ');
    writematrix(bias_lstm1,'bias_lstm1.txt','Delimiter',' ');
    writematrix(bias_lstm2,'bias_lstm2.txt','Delimiter',' ');
    writematrix(weights_cl,'weights_cl.txt','Delimiter',' ');
    writematrix(bias_cl,'bias_cl.txt','Delimiter',' ');
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
    networks = sprintf('network_%d_%d',length(cstate1),n);
    mkdir(networks);
    cd (networks)
    writematrix(muX,'muX.txt','Delimiter',' ');
    writematrix(muX,'muT.txt','Delimiter',' ');
    writematrix(sigmaX,'sigmaX.txt','Delimiter',' ');
    writematrix(sigmaT,'sigmaT.txt','Delimiter',' ');
    writematrix(iweights_lstm,'iweights_lstm.txt','Delimiter',' ');
    writematrix(rweights_lstm,'rweights_lstm.txt','Delimiter',' ');
    writematrix(bias_lstm,'bias_lstm.txt','Delimiter',' ');
    writematrix(weights_cl,'weights_cl.txt','Delimiter',' ');
    writematrix(bias_cl,'bias_cl.txt','Delimiter',' '); 
end

cd ../


end