function [Y, C, G] = lstmForward(X, W, R, b, c0, y0)
% lstmForward   Propagate Long Short Term Memory layer forwards on the host
%   [Y, C, G] = lstmForward(X, W, R, b, c0, y0) computes the forward
%   propagation of the Long Short Term Memory layer using input data X,
%   input weights W, recurrent weights R, bias term b and initial cell
%   state c0, and initial hidden units y0.
%
%   Definitions:
%   D := Number of dimensions of the input data
%   N := Number of input observations (mini-batch size)
%   S := Sequence length
%   H := Hidden units size
%
%   Inputs:
%   X - Input data            (D)x(N)x(S) array
%   W - Input weights         (4*H)x(D) matrix
%   R - Recurrent weights     (4*H)x(H) matrix
%   b - Bias                  (4*H)x(1) vector
%   c0 - Initial cell state   (H)x(1) vector
%   y0 - Initial hidden units (H)x(1) vector
%
%   Outputs:
%   Y - Output                (H)x(N)x(S) array
%   C - Cell state            (H)x(N)x(S) array
%   G - Gates                 (4*H)x(N)x(S) array

% Determine dimensions
%[~, N, S] = size(X);
N = 1;
S = 1;
H = size(R, 2);

% Pre-allocate output, gate vectors and cell state
G = zeros(4*H, N, S, 'like', X);
Y = zeros(H, N, S, 'like', X);
C = zeros(H, N, S, 'like', X);

% Indexing helpers
[zInd, iInd, fInd, oInd] = gateIndices(H);
ifoInd = [iInd fInd oInd];

% Forward propagate through time
% for tt = 1:S
%     if tt == 1
tt = 1;
        % Linear gate operations
        G(:, :, tt) = W*X(:, :, tt) + R*y0 + b;
        
        % Nonlinear gate operations
        G(zInd, :, tt) = tanh( G(zInd, :, tt) );
        G(ifoInd, :, tt) = sigmoidForward( G(ifoInd, :, tt) );
%         G = iNonlinearActivations( G, zInd, ifoInd, tt );
        
        % Cell state update
        C(:, :, tt) = G(zInd, :, tt) .*  G(iInd, :, tt) + ...
            G(fInd, :, tt) .* c0;
       
%     else
%         % Linear gate operations
%         G(:, :, tt) = W*X(:, :, tt) + R*Y(:, :, tt - 1) + b;
%         
%         % Nonlinear gate operations
%         G = iNonlinearActivations( G, zInd, ifoInd, tt );
%         
%         % Cell state update
%         C(:, :, tt) = G(zInd, :, tt) .*  G(iInd, :, tt) + ...
%             G(fInd, :, tt) .* C(:, :, tt - 1);
%         
%     end
    
    % Layer output
    Y(:, :, tt) = tanh( C(:, :, tt) ) .* G(oInd, :, tt);
        
end

function G = iNonlinearActivations( G, zInd, ifoInd, tt )
% Nonlinear gate operations
G(zInd, :, tt) = tanh( G(zInd, :, tt) );
G(ifoInd, :, tt) = sigmoidForward( G(ifoInd, :, tt) );
end

% function Z = tanhForward(X)
% % tanhForward   Tanh activation
% %
% % Input:
% % X - The input feature maps for a set of images. A (H)x(W)x(C)x(N) array.
% %
% % Output:
% % Z - The output feature maps for a set of images. A (H)x(W)x(C)x(N) array.
% 
% %   Copyright 2015-2016 The MathWorks, Inc.
% 
% Z = tanh(X);
% end

function Z = sigmoidForward(X)
% sigmoidForward   Sigmoid activation
%
% Input:
% X - The input feature maps for a set of images. A (H)x(W)x(C)x(N) array.
%
% Output:
% Z - The output feature maps for a set of images. A (H)x(W)x(C)x(N) array.

%   Copyright 2015-2016 The MathWorks, Inc.

Z = 1 ./ (1 + exp(-X));
end

function [zInd, iInd, fInd, oInd] = gateIndices(HiddenSize)
% gateIndices   Determine indices of the data input, input, forget and
% output gates of the LSTM layer

%   Copyright 2017 The MathWorks, Inc.

iInd = 1:HiddenSize;
fInd = 1 + HiddenSize:2*HiddenSize;
zInd = 1 + 2*HiddenSize:3*HiddenSize;
oInd = 1 + 3*HiddenSize:4*HiddenSize;

end