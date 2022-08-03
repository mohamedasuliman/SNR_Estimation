% Coded by: Mohamed A. Suliman.
% E-mail: mohamedabdall78@hotmail.com
% Date: April 15, 2017
% snr_estimation version 1.0
%   ----------------
%  | snr_estimation | is a Matlab function for estimating the SNR by estimating the in the linear system.
%   ----------------
% 
% This function implements the algorithm described in the paper: "SNR Estimation in Linear Systems", 
% Mohamed Suliman, Ayed M. Alrashdi, Tarig Ballal, & Tareq Y. Al-Naffouri,
% IEEE Signal Processing Letters, 24 (12), 1867-1871, April 2017.
%
%                            Brief Summary 
%                            -------------
% We estimate the SNR for the linear systems in the form
%                           y = Wx + n,
% Basic assumtions:
% -----------------
%
% 1) W  =  Phsi^(1/2)*W_bar.  (M by K)
% 2) W_bar : Gaussian matrix with i.i.d. entries of zero mean and unit variance.
% 3) The entries of x, and those of the noise n, are drawn from any two distributions 
%    and they are i.i.d. with zero mean and unknown variances.

%                              Function details 
%                              ----------------
% Inputs:
% ------- 
% y    :  Observations vector of length M.
% W_bar:  (M by K) Gaussian channel matrix with i.i.d. entries.
% Phsi :  (M by M) known Hermitian nonnegative left correlation matrix for W_bar.
% M    :  Number of observations.
% K    :  Number of unkowns.

% Outputs:
% --------

% SNR       :   The signal to noise ratio (SNR) of the linear system. 
% signal_var:   The signal varince.
% noise_var :   The noise varince.
% 

function  [SNR ,signal_var, noise_var]  = snr_estimation(y,W_bar,Psi,M,K)

lambda_values = .001:.001:.015;  % The values at which the deterministic equivalent will be evaluated
Num_lambdas   = numel(lambda_values);

I_K  = eye(K);
I_M  = eye(M);
[U,Q]= eig(Psi);
W    = Psi^(1/2)*W_bar;
b    = W'*y;
WW   = W'*W;

Xi_1 = zeros(Num_lambdas,1);
Xi_2 = zeros(Num_lambdas,1);
Phi  = zeros(Num_lambdas,1);

  for lambda_index = 1:1:Num_lambdas
    
     lambda = lambda_values(lambda_index);

     t= K/lambda;

    x_estimate  = (WW+lambda*I_K)^(-1)*b;  % RLS estimation for x.

    % Start solving the fixed point equation: delta (t) (Equation (11)).
    
    delta     =1;
    delta_aux =2;
    it        =0;

        while abs((delta-delta_aux)/delta)>10^-4
            it=it+1;
            delta_aux=delta;
            delta=1/K*trace(Q*inv(I_M+t/(1+t*delta_aux)*Q));      
        end

     T    = inv(I_M+t/(1+t*delta)*Psi);   % T(t) matrix (Equation (10))
     Temp = trace(Psi*T) ;                

     Xi_1(lambda_index,1) = (Temp/((1+t*delta))); 
     Xi_2(lambda_index,1) = (M/K-t*Temp/(K*(1+t*delta)));
     Phi (lambda_index,1) = (1/K)*(norm(y-W*x_estimate)^(2)+lambda*norm(x_estimate)^(2));


  end
 
 % Now solve the linear system: Xi*sigmas + epsilon = Phi  (Equation (14))

 Xi        = [Xi_1   Xi_2]; 
 sigmas    = Xi\Phi; 
 signal_var= sigmas(1);
 noise_var = sigmas(2);
 SNR       = signal_var/noise_var;
 
 
end
 
 
