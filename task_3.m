clc;

% Original message


msg = repmat([1 0 1 1 0 1], 1, 10);
Es_N0 = 2;

% CRC encode

crcCfg = crcConfig;
codeword = crcGenerate(msg.', crcCfg).';

disp("CRC encoded");
disp(codeword);

% Convolutional encode (Task 2 function, or matlab built-in function)
% compare with built in matlab function
k = 3;
trellis = poly2trellis(k, [7 5]);

convencoder = comm.ConvolutionalEncoder(trellis);


encoded_msg_Comm = convencoder([codeword, zeros(1,k-1)].').'; % padding with zeros


%disp("Convolutional encoded");
%disp(encoded_msg_Comm);


% Simulate channel — introduce errors

bpskmodulator = comm.BPSKModulator;
modulated_signal = bpskmodulator(encoded_msg_Comm.').';
awgnchan = comm.AWGNChannel("NoiseMethod","Signal to noise ratio (Es/No)","EsNo",Es_N0,"BitsPerSymbol",1);


noisy_signal = awgnchan(modulated_signal.').';

bpskdemodulator = comm.BPSKDemodulator;
recieved_signal = bpskdemodulator(noisy_signal.').';

% Viterbi decode

%viterbiDecoder = comm.ViterbiDecoder(TrellisStructure=trellis, InputFormat="Hard");
%decoded_msg = viterbiDecoder(noisy_signal.').';
decoded_msg = vitdec(recieved_signal, trellis, 2, "term", "hard");
decoded_msg = decoded_msg(:, 1 : end-(k-1));

disp("Viterbi decoded:");
disp(decoded_msg);

% CRC check — did Viterbi recover everything correctly?
[tx, err] = crcDetect(decoded_msg.', crcCfg);

disp("tx");
disp(tx.');
disp("err");
disp(err);

