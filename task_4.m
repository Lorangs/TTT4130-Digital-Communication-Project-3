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
iter_num = 200;
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
    %constraintLength = 7; 
    
    % 1. Return to the Rate 1/2 Mother Code
    %generatorPolys = [171 133]; 
    
    constraintLength = 9;
    generatorPolys = [561 753];
    
    message = randi([0,1],mess_length,1);
    poly = 'z^32 + z^26 + z^23 + z^22 + z^16 + z^12 + z^11 + z^10 + z^8 + z^7 + z^5 + z^4 + z^2 + z + 1';
    
    % Create the CRC Generator
    crcConf = crcConfig(Polynomial=poly, InitialConditions=1, DirectMethod=true,FinalXOR=1);
    MessageCrc = crcGenerate(message, crcConf);
    
    trellis = poly2trellis(constraintLength, generatorPolys);
    
    % 2. Puncture Pattern for Rate 7/10
    % A Rate 1/2 encoder turns 7 input bits into 14 output bits.
    % We want 10 output bits, so we delete (0) 4 bits out of every 14.
    % Distributing the zeros helps the decoder perform better.
    puncturePattern = [1; 1; 1; 1; 1; 1; 1; 0; 1; 1; 1; 1; 1; 0]; 
    
    encodedData = convenc(MessageCrc, trellis, puncturePattern);
    
    % Calculate the actual code rate to verify
    actualRate = length(MessageCrc) / length(encodedData);
    fprintf('Target Rate: 0.7000 (7/10)\n');
    fprintf('Actual Rate achieved: %f\n', actualRate);
    fprintf('Info Rate: %f bits/sym\n\n', actualRate * 4);
    
    scrambled_bits = matintrlv(encodedData, length(encodedData)/100,100);
    
    syms = qammod(scrambled_bits, 16, "gray", InputType="bit");
    
    % Channel
    snr_db = 10;
    RecSyms = awgn(syms, snr_db, "measured");
    
    % 3. SOFT DECISION DEMODULATION (Crucial for 10^-3 BER target)
    
    noiseVar = (10^(-snr_db/10));
    
    ReceivedData = qamdemod(RecSyms, 16, "gray", ...
        OutputType="approxllr", NoiseVariance=noiseVar);
    
    ReceivedData = matdeintrlv(ReceivedData, length(encodedData)/100, 100);
    
    % 4. SOFT DECISION VITERBI DECODING
    tracebackDepth = 7*constraintLength;
    
    % 'unquant' feeds the exact LLR math into the Viterbi decoder
    decodedData = vitdec(ReceivedData, trellis, tracebackDepth, 'trunc', 'unquant', puncturePattern);
    
    [rxDataBits, err] = crcDetect(decodedData, crcConf);
    
    % Calculate the Bit Error Rate (BER)
    [numErrors, BER] = biterr(message, rxDataBits(1:mess_length));
    errors_after_dec = numErrors;
    errors_crc = err;
end