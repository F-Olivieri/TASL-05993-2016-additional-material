%% Supplemental material for T-ASL-05993-2016
% Copyright (c) 2017, Ferdinando Olivieri,
% Institute of Sound and Vibration Research, University of Southampton, United Kingdom.
% E-mail: f.olivieri@ieee.org
% Last Update: 31st January 2017

% Permission to use, copy, modify, and/or distribute this software for any purpose
% with or without fee is hereby granted, provided that the above copyright notice
% and this permission notice appear in all copies.

% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL,
% DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
% RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
% NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE
% USE OR PERFORMANCE OF THIS SOFTWARE.

% A README file should have been provided with this code. If that is not the case,
% please send an email to f.olivieri@ieee.org and request one.
% The README file contains detailed information about this Matlab script and the associated data.

%% Description
% This data set and code accompanies the manuscript:
% Ferdinando Olivieri, Filippo Maria Fazi, Simone Fontana, Dylan Menzies and Philip Arthur Nelson,
% "Generation of private sound with a circular loudspeaker array and the Weighted Pressure Matching method ",
% submitted for consideration to IEEE Transactions on Audio, Speech and Language Processing
% (Reference Number T-ASL-05993-2016).
% DOI: <TBA>

% The code contains the implementations of Full-search of psi_D, Alg 2, and Alg 4 for the QCS scenario.
% Please refer to the information contained in the paper for a detailed explanation of each algorithm.

% The provided code has been tested with Matlab 2016b.
% The functions that implement the algorithm are provided in the FUNCTIONS
% section of this script

% When using any part of this code or data, please cite above paper.
% We would appreciate if you send us a note to the email address provided above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables initialization
clear vars; close all; clc; % Reset the workspace
load('OlivieriEtAl_TASL059932016_data.mat'); % The MAT file that contains the transfer functions Z, N, and the vector of the temporal frequency

controlpoint_angles = 0:5:355; % in degrees
freq_limits = [100, 10000]; % Frequencies in Hz. Limits of the x-asis
M = size(Z, 1); % Number of control points
L = size(Z, 2); % Number of loudspeakers
lengthfreqbins =  size(Z, 3); % Number of frequency bins

% The transfer functions at the control points
z_B = squeeze(Z(bp_idx, :, :));
Z_D = Z(dp_idx, :, :);

beta0 = 10^-2; % The value of $\beta_0$ (Tichonov Regularization Factor)
pB_min = 0.7079; % the -3dB condition (set by the user)

NUM_REPETITION = 100; % Number of repetitions for each algorithm (for
% the calculation of the average running time)

%% SYSTEM INITIALIZATION

% Calculate Regularization parameter by means of the the Normalized Tikhonov
% regularization technique (see Section III.B "Choice of the regularization
% parameter".
betavect = zeros(lengthfreqbins, 1);
for freq_idx = 1:lengthfreqbins
    Sigma = svd(Z(:, :, freq_idx));
    sigmasquared = max(Sigma(:))^2;
    betavect(freq_idx) = beta0*sigmasquared;
end

% Plots
figure; set(gcf,'color','w');
subplot(2,1,1);
semilogx(freqvect, betavect);
xlim(freq_limits); grid on; xlabel('Frequency, Hz');
ylabel('Regularization factor');

subplot(2,1,2);
semilogx(freqvect, N); ylim([0, 40]);
xlim(freq_limits); grid on; xlabel('Frequency, Hz'); ylabel('N');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Common Matrices for Full-search of psi_D and Alg 2
I = eye(L, L);
zBT= zeros(lengthfreqbins, L );
zBconj = zeros(L, lengthfreqbins);
zBconj_zBT =  zeros(L, L, lengthfreqbins);
ZDH_ZD =  zeros(L, L, lengthfreqbins);
betaI =  zeros(L, L, lengthfreqbins);
for freq_idx = 1:lengthfreqbins
    curr_ZB = z_B(:, freq_idx);
    curr_ZD = Z_D(: ,:, freq_idx);
    zBT(freq_idx, :) = transpose(curr_ZB);
    zBconj_zBT(:, :, freq_idx) = conj(curr_ZB)*transpose(curr_ZB);
    zBconj(:, freq_idx) = conj(curr_ZB);
    ZDH_ZD(:, :, freq_idx) =  ctranspose(curr_ZD)*curr_ZD;
    betaI(:, :, freq_idx) = betavect(freq_idx)*I;
end

%% 2) Extra matrices for Alg 2. Initialization (offline)
psi_D_hat = 0.5;
Aminus1 = zeros(L, L, lengthfreqbins);
Aminus1_ZHDZH = zeros(L, L, lengthfreqbins);
q_psiD_hat = zeros(L, lengthfreqbins); % The initial set of input signals Eq 18
for freq_idx = 2:lengthfreqbins
    Aminus1(:, :, freq_idx) = inv(zBconj_zBT(:, :, freq_idx) + psi_D_hat*ZDH_ZD(:, :, freq_idx) + betaI(:, :, freq_idx));
    Aminus1_ZHDZH(:, :, freq_idx) =  Aminus1(:, :, freq_idx)*ZDH_ZD(:, :, freq_idx); % For the calculation of pBn in Eq 26
    q_psiD_hat(:, freq_idx) = Aminus1(:, :, freq_idx)*zBconj(:, freq_idx); % Eq 18
end

% Calculate pB,n Eq 26 and store them in the system
pBn = cell(lengthfreqbins, 1);
for freq_idx = 1:lengthfreqbins
    pBn_currFreq = zeros(N(freq_idx), 1);
    idx = 1;
    for n = N(freq_idx):-1:0 % Store the elements of the polynomial from higher N to lower order
        pBn_currFreq(idx) = (-1)^n * zBT(freq_idx, :) * (Aminus1_ZHDZH(:, :, freq_idx))^n * q_psiD_hat(:, freq_idx);
        idx = idx + 1;
    end
    pBn{freq_idx} = pBn_currFreq;
end


%% Running Full-search of psi_D (FSP)
MatTemp = zeros(L, L, lengthfreqbins);
for freq_idx = 2:lengthfreqbins
    MatTemp(:, :, freq_idx) =  zBconj_zBT(:, :, freq_idx)  + betaI(:, :, freq_idx);
end

delta_psi_fsp = 10^-2; % Value used (see Experimental Validation Section)

q_fsp = zeros(L, lengthfreqbins);
psi_D_fsp = zeros(lengthfreqbins, 1);
% The full-search algorithm (implementation of PseudoCode 1)
disp('Running Full-search of psi_D... Please wait (it may take a long time...)');
tic;
[q_fsp, psi_D_fsp] = fsp(q_fsp, psi_D_fsp, MatTemp, ZDH_ZD, zBconj, zBT, pB_min, delta_psi_fsp, lengthfreqbins);
toc;
disp('Running Full-search of psi_D... Done!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running Alg 1
% Calculation (once user changes pB_min)
ImagPartThreshols = 10^-7;
PSI_MAX = 1; PSI_MIN = 0;
psi_D_Alg1 = zeros(lengthfreqbins, 1);
q_tilde_Alg1 = zeros(L, lengthfreqbins);
t_alg1 = zeros(1, NUM_REPETITION);
disp(['Running Alg 2 for ' num2str(NUM_REPETITION) ' iterations... Please wait']);
for n = 1:NUM_REPETITION
    tic;
    [q_tilde_Alg1, psi_D_Alg1] = Alg1(q_tilde_Alg1, psi_D_Alg1, MatTemp, ...
        ZDH_ZD, zBconj, pBn, N, pB_min, psi_D_hat, PSI_MIN, PSI_MAX, ImagPartThreshols, lengthfreqbins);
    t_alg1(n) = toc;
end
% toc;
disp('... Done!');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running Alg 2

Psi_Max_INIT = 1;
Psi_Min_INIT = 0;
thresh = 10^-7;

q_Alg2 = zeros(L, lengthfreqbins);
psi_D_Alg2 = zeros(lengthfreqbins, 1);
t_alg2 = zeros(1, NUM_REPETITION);
disp(['Running Alg 4 for ' num2str(NUM_REPETITION) ' iterations... Please wait']);
for n = 1:NUM_REPETITION
    tic;
    for freq_idx = 2:lengthfreqbins
        [psi_D_Alg2(freq_idx), q_Alg2(:, freq_idx)] = Alg2_RecursiveBisectionSearch(zBT(freq_idx, :), zBconj_zBT(:, :, freq_idx), ...
            betaI(:, :, freq_idx), ZDH_ZD(:, :, freq_idx), zBconj(:, freq_idx), pB_min, ...
            Psi_Max_INIT, Psi_Min_INIT, thresh);
    end
    t_alg2(n) = toc;
end
disp('... Done!');

%% Comparison between running times (Alg 1 and Alg 2)

figure; set(gcf,'color','w');
plot(t_alg1, '-')
hold all; plot(t_alg2, ':'); grid on;
legend('Alg 2', 'Alg 4', 'Location', 'Best');
xlabel('Repetition'); ylabel('Time in seconds');
xlim([1, 1000]); ylim([0, 4]);
disp(['Average time Alg 2:' num2str(mean(t_alg1)) ' s']);
disp(['Average time Alg 4:' num2str(mean(t_alg2)) ' s']);

%% Plot of the psi_D

figure;
semilogx(freqvect, psi_D_fsp);
hold all; semilogx(freqvect, psi_D_Alg1);
hold all; semilogx(freqvect, psi_D_Alg2);
legend('Full-search', 'Alg 1', 'Alg 2');
xlim(freq_limits); ylim([0, 0.4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterative search
function [q_fsp, psi_D_fsp] = fsp(q_fsp, psi_D_fsp, MatTemp, ZDH_ZD, zBconj, zBT, pB_min, delta_psi_Alg_1, lengthfreqbins, L)
% Implementation of Full-search of psi_D
% Author: F Olivieri

for freq_idx = 2:lengthfreqbins
    psi_D_temp = 1;
    q_temp = (MatTemp(:, :, freq_idx) + psi_D_temp*ZDH_ZD(:, :, freq_idx))\ zBconj(:, freq_idx);
    
    while zBT(freq_idx, :)*q_temp < pB_min && psi_D_temp > 0
        psi_D_temp = psi_D_temp - delta_psi_Alg_1;
        q_temp = (MatTemp(:, :, freq_idx) + psi_D_temp*ZDH_ZD(:, :, freq_idx))\ zBconj(:, freq_idx);
    end
    q_fsp(:, freq_idx) = q_temp;
    psi_D_fsp(freq_idx) = psi_D_temp;
end

end

%% Algorithm 1
function [q_tilde_Alg1,psi_D_Alg1 ] = Alg1(q_tilde_Alg1, psi_D_Alg1, MatTemp, ...
    ZDH_ZD, zBconj, pBn, N, pB_min, psi_D_hat, PSI_MIN, PSI_MAX, ImagPartThreshols, lengthfreqbins)
% Implementation of Alg 2
% Author: F Olivieri
for freq_idx = 2:lengthfreqbins
    polyB = pBn{freq_idx};
    
    last_elem_idx = N(freq_idx) + 1;
    
    % Eq 28
    polyB(last_elem_idx) = polyB(last_elem_idx) - pB_min; % subtract pB_min from the term of order zero
    rootsB_ = roots(polyB); % Calculating the roots of the polynomial
    est_AbsDeltaPsi0_B = real(rootsB_(abs(imag(rootsB_)) < ImagPartThreshols)); % keep roots whose imaginary parts is approx zero
    est_Psi0 = psi_D_hat + est_AbsDeltaPsi0_B;
    
    % Eq 29
    if est_Psi0 < PSI_MIN, est_Psi0 = PSI_MIN;
    elseif est_Psi0 > PSI_MAX, est_Psi0 = PSI_MAX;
    end
    
    psi_D_Alg1(freq_idx) = est_Psi0;
    
    q_tilde_Alg1(:, freq_idx) = (MatTemp(:, :, freq_idx) + est_Psi0*ZDH_ZD(:, :, freq_idx))\ zBconj(:, freq_idx);
end

end

%% Algorithm 2
function [psi_est, q_temp] = Alg2_RecursiveBisectionSearch(zBT, zBT_zBconj, betaI, ZDH_ZD, zBconj, pBmin, Psi_Max, Psi_Min, thresh)
% Implementation of Alg 4
% Author: F Olivieri

psi_est = (Psi_Max - Psi_Min)/2 + Psi_Min;

q_temp = (zBT_zBconj + psi_est*ZDH_ZD + betaI)\zBconj;
pB = real(zBT*q_temp);
condition1 = abs(pB - pBmin);

if condition1 <= thresh ||  Psi_Max <= thresh % condition is met: Stop!
    return;
else % search for the parameter
    if pB < pBmin % yes: need to decrease psi_est
        Psi_Max =  psi_est;
    else % no: need to increase psi_est
        Psi_Min = psi_est;
    end
    
    psi_est = Alg2_RecursiveBisectionSearch(zBT, zBT_zBconj, betaI, ZDH_ZD, ...
        zBconj, pBmin, Psi_Max, Psi_Min, thresh);
end

end