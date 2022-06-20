function [Zica, W, T, mu, Zcw, normRows] = fastICA(Z,r,type,flag)
%
% Syntax:       Zica = fastICA(Z,r);
%               Zica = fastICA(Z,r,type);
%               Zica = fastICA(Z,r,type,flag);
%               [Zica, W, T, mu] = fastICA(Z,r);
%               [Zica, W, T, mu] = fastICA(Z,r,type);
%               [Zica, W, T, mu] = fastICA(Z,r,type,flag);
%               
% Inputs:       Z is an d x n matrix containing n samples of d-dimensional
%               data
%               
%               r is the number of independent components to compute
%               
%               [OPTIONAL] type = {'kurtosis','negentropy'} specifies
%               which flavor of non-Gaussianity to maximize. The default
%               value is type = 'kurtosis'
%               
%               [OPTIONAL] flag determines what status updates to print
%               to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      Zica is an r x n matrix containing the r independent
%               components - scaled to variance 1 - of the input samples
%               
%               W and T are the ICA transformation matrices such that
%               Zr = T \ W' * Zica + repmat(mu,1,n);
%               is the r-dimensional ICA approximation of Z
%               
%               mu is the d x 1 sample mean of Z
%               
% Description:  Performs independent component analysis (ICA) on the input
%               data using the Fast ICA algorithm
%               
% Reference:    Hyv?inen, Aapo, and Erkki Oja. "Independent component
%               analysis: algorithms and applications." Neural networks
%               13.4 (2000): 411-430
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 26, 2015
%               November 12, 2016
%

% Constants
TOL = 1e-6;         % Convergence criteria
MAX_ITERS = 100;    % Max # iterations

% Parse inputs
if ~exist('flag','var') || isempty(flag)
    % Default display flag
    flag = 1;
end
if ~exist('type','var') || isempty(type)
    % Default type
    type = 'kurtosis';
end

% Set algorithm type
if strncmpi(type,'kurtosis',1)
    % Kurtosis
    USE_KURTOSIS = true;
    algoStr = 'kurtosis';
elseif strncmpi(type,'negentropy',1)
    % Negentropy
    USE_KURTOSIS = false;
    algoStr = 'negentropy';
else
    % Unsupported type
    error('Unsupported type ''%s''',type);
end

% Center and whiten data
[Zc, mu] = centerRows(Z);
[Zcw, T] = whitenRows(Zc);

% Normalize rows to unit norm
normRows = @(X) bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));

% Perform Fast ICA
if flag
    % Prepare status updates
    fmt = sprintf('%%0%dd',ceil(log10(MAX_ITERS + 1)));
    str = sprintf('Iter %s: max(1 - |<w%s, w%s>|) = %%.4g\\n',fmt,fmt,fmt);
    fprintf('***** Fast ICA (%s) *****\n',algoStr);
end
W = normRows(rand(r,size(Z,1))); % Random initial weights
% W = normRows(eye(r,size(Z,1))); % Random initial weights
k = 0;
delta = inf;
while delta > TOL && k < MAX_ITERS
    k = k + 1;
    
    % Update weights
    Wlast = W; % Save last weights
    Sk = permute(W * Zcw,[1, 3, 2]);
    if USE_KURTOSIS
        % Kurtosis
        G = 4 * Sk.^3;
        Gp = 12 * Sk.^2;
    else
        % Negentropy
        G = Sk .* exp(-0.5 * Sk.^2);
        Gp = (1 - Sk.^2) .* exp(-0.5 * Sk.^2);
    end
    W = mean(bsxfun(@times,G,permute(Zcw,[3, 1, 2])),3) - ...
             bsxfun(@times,mean(Gp,3),W);

    W = normRows(W);
    
    % Decorrelate weights
    [U, S, ~] = svd(W,'econ');
    W = U * diag(1 ./ diag(S)) * U' * W;
    
    % Update convergence criteria
    delta = max(1 - abs(dot(W,Wlast,2)));
    if flag
        fprintf(str,k,k,k - 1,delta);
    end
end
if flag
    fprintf('\n');
end

% Independent components
Zica = W * Zcw;
end
%%
function [Zc, mu] = centerRows(Z)
%
% Syntax:       [Zc, mu] = centerRows(Z);
%               
% Inputs:       Z is an (d x n) matrix containing n samples of a
%               d-dimensional random vector
%               
% Outputs:      Zc is the centered version of Z
%               
%               mu is the (d x 1) sample mean of Z
%               
% Description:  Returns the centered (zero mean) version of the input data
%               
% Note:         Z = Zc + repmat(mu,1,n)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         November 1, 2016
%

% Compute sample mean
mu = mean(Z,2);

% Subtract mean
Zc = bsxfun(@minus,Z,mu);
end
%%
function [Zw, T] = whitenRows(Z)
%
% Syntax:       [Zw, T] = whitenRows(Z);
%               
% Inputs:       Z is an (d x n) matrix containing n samples of a
%               d-dimensional random vector
%               
% Outputs:      Zw is the whitened version of Z
%               
%               T is the (d x d) whitening transformation of Z
%               
% Description:  Returns the whitened (identity covariance) version of the
%               input data
%               
% Notes:        (a) Must have n >= d to fully whitenRows Z
%               
%               (b) Z = T \ Zcw
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         November 1, 2016
%

% Compute sample covariance
R = cov(Z');

% Whiten data
[U, S, V] = svd(R);
T  = U * diag(1 ./ sqrt(diag(S))) * U';
Zw = T * Z;
end