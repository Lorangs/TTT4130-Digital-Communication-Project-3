% Mutual Information Calculation for Discrete Channels (BPSK/QAM)
% with Hard Decision Decoding (Quantized Output)
clear; close all; clc;

%% 1. Configuration
bits_per_symbol_list = 2:2:8; % Order: 1=BPSK, 2=QPSK, 3=8QAM, 4=16QAM
SNR_dB_range = -10:0.5:30;  % SNR range in dB

% Initialize storage for plotting
MI_results = zeros(length(bits_per_symbol_list), length(SNR_dB_range));

%% 2. Main Loop
fprintf('Starting calculation...\n');

for i = 1:length(bits_per_symbol_list)
    k = bits_per_symbol_list(i);
    M = 2^k; % Constellation size
    
    % --- Generate Constellation X ---
    if k == 1
        % BPSK (Real only)
        X = [-1, 1];
        is_complex = false;
    else
        % QAM (Complex)
        is_complex = true;
        if k == 3
            % Rectangular 8-QAM (4x2 grid) to keep regions rectangular
            % This makes the analytical integral separable and exact.
            x_vals = [-(3):2:3]; % -3, -1, 1, 3
            y_vals = [-1, 1];    % -1, 1
            [XI, XQ] = meshgrid(x_vals, y_vals);
            X = XI(:)' + 1j*XQ(:)';
        else
            % Square QAM (4-QAM, 16-QAM)
            X = qammod(0:M-1, M, 'UnitAveragePower', false);
            % Ensure it is centered and square scaling is standard
            % (Standard qammod output is suitable for decision boundary calc)
        end
    end
    
    % Calculate Average Symbol Energy (Es)
    Es = mean(abs(X).^2);
    
    % --- Determine Decision Boundaries ---
    % Since Y is the discretization of T to the same values of X (Minimum Distance),
    % the boundaries are the midpoints between sorted unique coordinate values.
    
    if is_complex
        real_pts = unique(real(X));
        imag_pts = unique(imag(X));
        bounds_I = [-inf, (real_pts(1:end-1) + real_pts(2:end))/2, inf];
        bounds_Q = [-inf, (imag_pts(1:end-1) + imag_pts(2:end))/2, inf];
    else
        real_pts = unique(X);
        bounds_I = [-inf, (real_pts(1:end-1) + real_pts(2:end))/2', inf];
        bounds_Q = []; % Not used for BPSK
    end
    
    % --- SNR Loop ---
    for j = 1:length(SNR_dB_range)
        snr_db = SNR_dB_range(j);
        snr_lin = 10^(snr_db/10);
        
        % Calculate Noise Sigma
        % SNR = Es / sigma_total^2
        sigma_total = sqrt(Es / snr_lin);
        
        if is_complex
            % For complex, noise power splits between I and Q
            sigma_dim = sigma_total / sqrt(2);
        else
            % For real (BPSK), all noise is in the real dimension
            sigma_dim = sigma_total;
        end
        
        % --- Calculate Transition Matrix P(y|x) ---
        % P_yx(row, col) = P(Y = y_row | X = x_col)
        P_yx = zeros(M, M);
        
        for tx_idx = 1:M
            x_curr = X(tx_idx);
            
            for rx_idx = 1:M
                y_target = X(rx_idx);
                
                % Find integration limits for this target symbol Y
                % 1. Real Dimension Limits
                [~, idx_I] = min(abs(real_pts - real(y_target)));
                lim_I_low = bounds_I(idx_I);
                lim_I_high = bounds_I(idx_I+1);
                
                % Probability in I (using normcdf)
                % P = Phi((Upper-Mean)/sig) - Phi((Lower-Mean)/sig)
                prob_I = normcdf(lim_I_high, real(x_curr), sigma_dim) - ...
                         normcdf(lim_I_low, real(x_curr), sigma_dim);
                     
                % 2. Imaginary Dimension Limits (if complex)
                if is_complex
                    [~, idx_Q] = min(abs(imag_pts - imag(y_target)));
                    lim_Q_low = bounds_Q(idx_Q);
                    lim_Q_high = bounds_Q(idx_Q+1);
                    
                    prob_Q = normcdf(lim_Q_high, imag(x_curr), sigma_dim) - ...
                             normcdf(lim_Q_low, imag(x_curr), sigma_dim);
                else
                    prob_Q = 1; % No imaginary component
                end
                
                % Total Probability
                P_yx(rx_idx, tx_idx) = prob_I * prob_Q;
            end
        end
        
        % --- Calculate Mutual Information I(X;Y) ---
        % Assume equiprobable input P(x) = 1/M
        P_x = 1/M;
        
        % P(y) = sum over x of P(y|x)P(x)
        P_y = sum(P_yx .* P_x, 2); % Sum columns for each row
        
        % MI = sum_x sum_y P(x,y) log2( P(x,y) / (P(x)P(y)) )
        %    = sum_x sum_y P(y|x)P(x) log2( P(y|x) / P(y) )
        
        MI_val = 0;
        for tx = 1:M
            for rx = 1:M
                p_trans = P_yx(rx, tx);
                if p_trans > 1e-12 % Avoid log(0)
                    MI_val = MI_val + p_trans * P_x * log2(p_trans / P_y(rx));
                end
            end
        end
        
        MI_results(i, j) = MI_val;
    end
end

fprintf('Calculation complete.\n');

%% 3. Plotting
figure('Color', 'w');
styles = {'-', '--', '-.', '-'};

hold on;
for i = 1:length(bits_per_symbol_list)
    plot(SNR_dB_range, MI_results(i, :), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('%d-QAM/PSK (%d bits/sym)', 2^bits_per_symbol_list(i), bits_per_symbol_list(i)));
end

% Formatting
grid on;
xlabel('SNR [dB] = 10 log_{10}(E_s / \sigma_{total}^2)');
ylabel('Mutual Information [bit]');
title('Mutual Information vs SNR for Discrete Quantized Channels');
legend('Location', 'SouthEast');
ylim([0 4.5]);
xlim([min(SNR_dB_range) max(SNR_dB_range)]);
ax = gca();
ax.FontSize = 30;
legend(FontSize=30)

% Add Channel Capacity Reference (log2(1+SNR)) for comparison?
% Optional: Uncomment to see Shannon Capacity for continuous AWGN
 cap_awgn = log2(1 + 10.^(SNR_dB_range./10));
 plot(SNR_dB_range, cap_awgn, 'k:', 'LineWidth', 1.5, 'DisplayName', 'AWGN Capacity (log_2(1+SNR))');

hold off;