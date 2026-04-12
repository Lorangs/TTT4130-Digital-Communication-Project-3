
clearvars
close all

snr_db = 0:1:15;
M = 16;
k = log2(M);
num_words = 10000;
num_bits = num_words*k;
msg = randi([0 1], num_bits, 1);

%% Uncoded message
uncoded_msg = msg;

numErrors_vec_uncoded = zeros(1, length(snr_db));
ber_vec_uncoded = zeros(1, length(snr_db));
modulated_uncoded_msg = qammod(uncoded_msg, M, "gray","InputType","bit");
for i = 1:length(snr_db)

    noisy_uncoded_msg = awgn(modulated_uncoded_msg, snr_db(i), "measured");
    demodulated_uncoded_msg = qamdemod(noisy_uncoded_msg, M, "gray", "OutputType","bit");

    % Compute BER for this SNR
    [numErrors_vec_uncoded(1,i), ber_vec_uncoded(1,i)] = biterr(msg(:), demodulated_uncoded_msg(:));

end

%% Hamming message
% Hamming code

n = 7;
m = n - k;

% Generate different seeds for random functions
rng('shuffle');


[H, G] = hammgen(m);
LUT = syndtable(H);
SyndTable = mod(LUT * H.', 2);

numErrors_vec_hamming = zeros(1, length(snr_db));
ber_vec_hamming = zeros(1, length(snr_db));
hamming_encoded_msg = zeros(n*num_words,1);
for j = 1:num_words
    hamming_encoded_msg((j-1)*n +1 : j*n) = encode(msg((j-1)*k +1 : j*k), n, k);
end
modulated_hamming_msg = qammod(hamming_encoded_msg, M, "gray","InputType","bit").';

for i = 1:length(snr_db)

    noisy_hamming_msg = awgn(modulated_hamming_msg, snr_db(i), "measured");
    demodulated_hamming_msg = qamdemod(noisy_hamming_msg, M, "gray", "OutputType","bit");
    
    decoded_hamming_msg = zeros(1, num_bits);
    for j = 1:num_words
        r = demodulated_hamming_msg((j-1)*n +1 : j*n);
        s = mod(r * H.', 2);
        [~, calc_e_idx] = ismember(s, SyndTable, 'rows');
        r = mod(r + LUT(calc_e_idx, :), 2);
        decoded_hamming_msg((j-1)*k +1 : j*k) = decode(r, n, k);
    end

    % Compute BER for this SNR
    [numErrors_vec_hamming(1,i), ber_vec_hamming(1,i)] = biterr(msg(:), decoded_hamming_msg(:));

end

%% Convolutional hard

k_conv = 7;
generatorPolys = [171 133]; 
trellis = poly2trellis(k_conv, generatorPolys);
convencoder = comm.ConvolutionalEncoder(trellis);

numErrors_vec_conv_hard = zeros(1, length(snr_db));
ber_vec_conv_hard = zeros(1, length(snr_db));
encoded_msg_conv_hard = convenc(msg, trellis); % padding with zeros
modulated_signal = qammod(encoded_msg_conv_hard, M, "gray","InputType","bit");
for i = 1:length(snr_db)
    
    noisy_signal = awgn(modulated_signal, snr_db(i), "measured");
    demodulated_msg_conv_hard = qamdemod(noisy_signal, M, "gray", "OutputType","bit");

    decoded_msg_conv_hard = vitdec(demodulated_msg_conv_hard, trellis, 5*k_conv, "trunc", "hard");
    %decoded_msg_conv_hard = decoded_msg_conv_hard(1 : end-(k_conv-1)); % remove padding

    % Compute BER for this SNR
    [numErrors_vec_conv_hard(i), ber_vec_conv_hard(i)] = biterr(msg(:), decoded_msg_conv_hard(:));
end

%% Convolutional soft

numErrors_vec_conv_soft = zeros(1, length(snr_db));
ber_vec_conv_soft = zeros(1, length(snr_db));

encoded_msg_conv_soft = convenc(msg, trellis); % padding with zeros
modulated_signal = qammod(encoded_msg_conv_soft, M, "gray","InputType","bit");
    
for i = 1:length(snr_db)
    
    noisy_signal = awgn(modulated_signal, snr_db(i), "measured");
    demodulated_msg_conv_soft = qamdemod(noisy_signal, M, "gray", "OutputType","approxllr", "NoiseVariance",10^(-snr_db(i)/10));

    decoded_msg_conv_soft = vitdec(demodulated_msg_conv_soft, trellis,5*k_conv, "trunc", 'unquant');
    %decoded_msg_conv_soft = decoded_msg_conv_soft(1 : end-(k_conv-1)); % remove padding

    % Compute BER for this SNR
    [numErrors_vec_conv_soft(1,i), ber_vec_conv_soft(1,i)] = biterr(msg(:), decoded_msg_conv_soft(:));
end

%%
figure;
semilogy(snr_db, ber_vec_hamming, '-o');
hold on;
semilogy(snr_db, ber_vec_uncoded, '-o');
semilogy(snr_db, ber_vec_conv_hard, '-o');
semilogy(snr_db, ber_vec_conv_soft, '-o');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for Hamming Code');
legend('Hamming', "Uncoded", 'Conv coder Hard', 'Conv coder soft');
grid on;
