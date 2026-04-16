clear all
close all
frame_length = 20*51;

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
fprintf('Total number of errors: %f ', sum(errors_after_dec));
fprintf('BER: %f\n', sum(errors_after_dec)/(mess_length*iter_num));

function [BER, errors_after_dec, errors_crc, mess_length] = coding_10db(frame_length)

    crc_length = 32;
    mess_length = frame_length - crc_length;
    constraintLength = 14; 
    
    % 1. Use a Rate 1/3 Mother Code to allow targeting ~1.7 bits/sym
    generatorPolys = [21645 35661 37133]; 
    
    message = randi([0,1],mess_length,1);
    poly = 'z^32 + z^26 + z^23 + z^22 + z^16 + z^12 + z^11 + z^10 + z^8 + z^7 + z^5 + z^4 + z^2 + z + 1';
    
    % Create the CRC Generator
    crcConf = crcConfig(Polynomial=poly, InitialConditions=1, DirectMethod=true,FinalXOR=1);
    MessageCrc = crcGenerate(message, crcConf);
    
    trellis = poly2trellis(constraintLength, generatorPolys);
   
    puncVector = ones(51, 1);
    puncIndices = 3:5:(5 * 5);
    puncVector(puncIndices) = 0;
    
    encodedData = convenc(MessageCrc, trellis, puncVector);
    actualRate = length(MessageCrc) / length(encodedData);
    fprintf('Target Rate: %f (17/40)\n', 17/40);
    fprintf('Actual Rate achieved: %f\n', actualRate);
    fprintf('Info Rate: %f bits/sym\n\n', actualRate * 4);
    
    sep_len = 120;
    scrambled_bits = matintrlv(encodedData, length(encodedData)/sep_len,sep_len);
    syms = qammod(scrambled_bits, 16, "gray", InputType="bit");
    
    % Channel
    snr_db = 5;
    RecSyms = awgn(syms, snr_db, "measured");
    
    % 3. SOFT DECISION DEMODULATION (Crucial for low BER)
    % Standard 16-QAM power is 10. We calculate exact noise variance for the LLRs. 
    noiseVar = (10^(-snr_db/10));
    
    % OutputType="approxllr" generates Soft Decisions
    ReceivedData = qamdemod(RecSyms, 16, "gray", ...
        OutputType="approxllr", NoiseVariance=noiseVar);
    
    % Deinterleave the soft LLR values
    ReceivedData = matdeintrlv(ReceivedData, length(encodedData)/sep_len, sep_len);
    
    % 4. SOFT DECISION VITERBI DECODING
    tracebackDepth = 7*constraintLength; % Keeps it a clean multiple of the puncture period (7)
    
    % Use 'unquant' to feed the LLRs into the Viterbi decoder
    decodedData = vitdec(ReceivedData, trellis, tracebackDepth, 'trunc', 'unquant', puncVector);
    
    [rxDataBits, err] = crcDetect(decodedData, crcConf);
    
    % Calculate the Bit Error Rate (BER)
    [numErrors, BER] = biterr(message, rxDataBits(1:mess_length));
    errors_after_dec = numErrors;
    errors_crc = err;

end