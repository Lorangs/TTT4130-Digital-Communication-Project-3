clc;

snr_db = 0:10;

k = 4;
num_words = 3000;
num_bits = num_words*k;
msg = randi([0 1], 1, num_bits);



%% Uncoded message
uncoded_msg = msg;

numErrors_vec_uncoded = zeros(1, length(snr_db));
ber_vec_uncoded = zeros(1, length(snr_db));

for i = 1:length(snr_db)
    awgnchan = comm.AWGNChannel("NoiseMethod","Signal to noise ratio (SNR)","SNR",snr_db(i),"BitsPerSymbol",1);


    modulated_uncoded_msg = pskmod(uncoded_msg.', 2).';
    noisy_uncoded_msg = awgnchan(modulated_uncoded_msg.').';
    demodulated_uncoded_msg = pskdemod(noisy_uncoded_msg.', 2).';

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

for i = 1:length(snr_db)
    awgnchan = comm.AWGNChannel("NoiseMethod","Signal to noise ratio (SNR)","SNR",snr_db(i),"BitsPerSymbol",1);
    hamming_encoded_msg = zeros(1, n*num_words);
    for j = 1:num_words
        hamming_encoded_msg((j-1)*n +1 : j*n) = encode(msg((j-1)*k +1 : j*k), n, k);
    end

    modulated_hamming_msg = pskmod(hamming_encoded_msg.', 2).';
    noisy_hamming_msg = awgnchan(modulated_hamming_msg.').';
    demodulated_hamming_msg = pskdemod(noisy_hamming_msg.', 2).';

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

k_conv = 3;
trellis = poly2trellis(k_conv, [7 5]);
convencoder = comm.ConvolutionalEncoder(trellis);

numErrors_vec_conv_hard = zeros(1, length(snr_db));
ber_vec_conv_hard = zeros(1, length(snr_db));

for i = 1:length(snr_db)
    awgnchan = comm.AWGNChannel("NoiseMethod","Signal to noise ratio (SNR)","SNR",snr_db(i),"BitsPerSymbol",1);
    encoded_msg_conv_hard = convencoder([msg, zeros(1,k_conv-1)].').'; % padding with zeros
    modulated_signal = pskmod(encoded_msg_conv_hard.', 2).';
    noisy_signal = awgnchan(modulated_signal.').';
    demodulated_msg_conv_hard = pskdemod(noisy_signal.', 2).';

    decoded_msg_conv_hard = vitdec(demodulated_msg_conv_hard, trellis, 2, "term", "hard");
    decoded_msg_conv_hard = decoded_msg_conv_hard(1, 1 : end-(k_conv-1)); % remove padding

    % Compute BER for this SNR
    [numErrors_vec_conv_hard(1,i), ber_vec_conv_hard(1,i)] = biterr(msg(:), decoded_msg_conv_hard(:));
end

%% Convolutional soft (not finished)

numErrors_vec_conv_soft = zeros(1, length(snr_db));
ber_vec_conv_soft = zeros(1, length(snr_db));

vitdec_object = comm.ViterbiDecoder(...
    TrellisStructure=trellis, ...
    InputFormat="Soft", ...
    TracebackDepth=2*k_conv,...
    TerminationMethod="Terminated" ...
);

for i = 1:length(snr_db)
    awgnchan = comm.AWGNChannel("NoiseMethod","Signal to noise ratio (SNR)","SNR",snr_db(i),"BitsPerSymbol",1);
    encoded_msg_conv_soft = convencoder([msg, zeros(1,k_conv-1)].').'; % padding with zeros
    modulated_signal = pskmod(encoded_msg_conv_soft.', 2).';
    noisy_signal = awgnchan(modulated_signal.').';
    demodulated_msg_conv_soft = pskdemod(noisy_signal.', 2, pi, "OutputType","approxllr").';

    decoded_msg_conv_soft = vitdec_object(demodulated_msg_conv_soft.').';
    decoded_msg_conv_soft = decoded_msg_conv_soft(1, 1 : end-(k_conv-1)); % remove padding

    % Compute BER for this SNR
    [numErrors_vec_conv_soft(1,i), ber_vec_conv_soft(1,i)] = biterr(msg(:), decoded_msg_conv_soft(:));
end






%%
figure;
semilogy(snr_db, ber_vec_hamming, '-o');
hold on;
semilogy(snr_db, ber_vec_uncoded, '-o');
hold on;
semilogy(snr_db, ber_vec_conv_hard, '-o');
hold on;
semilogy(snr_db, ber_vec_conv_soft, '-o');
hold off;
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER vs SNR for Hamming Code (Log-Log Scale)');
legend('Hamming', "Uncoded", 'Conv coder Hard', 'Conv coder soft');
grid on;



