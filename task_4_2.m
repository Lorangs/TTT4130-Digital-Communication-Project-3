clear all
close all
frame_length = 14*100;

%% 16QAM channel capacity 
MI_16QAM_10db = 2.8;
MI_16QAM_5db = 1.7;
data_rate_16QAM_10db = 2.8/4;
data_rate_16QAM_5db = 1.7/4;

MI_256QAM_10db = 3.2;
MI_256QAM_5db = 1.9;
data_rate_256QAM_10db = 3.2/8;
data_rate_256QAM_5db = 1.9/8;
BER = 0;
iter_num = 20000;
errors_after_dec = zeros(1, iter_num);
errors_crc = zeros(1, iter_num);

for i = 1:iter_num
    [BER, errors_after_dec(i), errors_crc(i), mess_length] = coding_10db(frame_length);
end
fprintf('Total number of errors: %f', sum(errors_after_dec));
fprintf('BER: %f\n', sum(errors_after_dec)/(mess_length*iter_num));

%% Function 10 dB SNR

function [BER, errors_after_dec, errors_crc, mess_length] = coding_10db(frame_length)

crc_length = 32;
mess_length = frame_length - crc_length;
constraintLength = 7; 
% Standard polynomials used in Wi-Fi/LTE (in octal)
generatorPolys = [171 133]; 
message = randi([0,1],mess_length,1);
poly = 'z^32 + z^26 + z^23 + z^22 + z^16 + z^12 + z^11 + z^10 + z^8 + z^7 + z^5 + z^4 + z^2 + z + 1';

% Create the CRC Generator
crcConf = crcConfig(Polynomial=poly, InitialConditions=1, DirectMethod=true,FinalXOR=1);
crcGen = crcGenerate(message, crcConf);
MessageCrc = crcGen;

trellis = poly2trellis(constraintLength, generatorPolys);

% 2. Define the Puncture Pattern for Rate 5/7
% We want 5 input bits to result in 7 output bits.
% A Rate 1/2 encoder turns 5 input bits into 10 output bits.
% Therefore, we need to delete (puncture) 3 out of every 10 bits.
% 1 = transmit, 0 = puncture (delete)
puncturePattern = [1; 1; 1; 0; 1; 1; 1; 0; 1; 0; 1; 1; 1; 0]; 

encodedData = convenc(MessageCrc, trellis, puncturePattern);

% Calculate the actual code rate to verify
actualRate = length(message) / length(encodedData);
fprintf('Target Rate: %f (5/7)\n', 5/7);
fprintf('Actual Rate achieved: %f\n', actualRate);

syms = qammod(encodedData, 4, "gray", InputType="bit", UnitAveragePower=true);

RecSyms = awgn(syms, 10);

ReceivedData = qamdemod(RecSyms, 4, "gray",OutputType="bit");

% 6. Decode the Data (Using Hard Decision Viterbi)
tracebackDepth = 32; % A standard rule of thumb is 5x constraint length
% We must pass the exact same puncture pattern to the decoder so it can
% insert "dummy" metrics where the bits were deleted.
decodedData = vitdec(ReceivedData, trellis, tracebackDepth, 'trunc', 'hard', puncturePattern);

[rxDataBits, err] = crcDetect(decodedData, crcConf);

% 7. Calculate the Bit Error Rate (BER)
% Compare the original data to the decoded data
[numErrors, BER] = biterr(message, rxDataBits(1:mess_length));
errors_after_dec = numErrors;
errors_crc = err;
end