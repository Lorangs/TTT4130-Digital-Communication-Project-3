clearvars
close all

SNR = 0:1:20;

message_length = 100;
constraintLength = 7; 
% Standard polynomials used in Wi-Fi/LTE (in octal)
generatorPolys = [171 133];

poly = 'z^32 + z^26 + z^23 + z^22 + z^16 + z^12 + z^11 + z^10 + z^8 + z^7 + z^5 + z^4 + z^2 + z + 1';
trellis = poly2trellis(constraintLength, generatorPolys);

k_ham = 10;
h = hammgen(5, [1,0,1,1]);

for i = 1:length(SNR)

    mess = randi([0,1],message_length,1);
    mess_ham = reshape(mess, [message_length/k_ham, k_ham]);

    hamm_encoded = encode(mess_ham, 5, k_ham, 'hamming/binary');
    conv_encoded = convenc(mess, trellis);

    sym_hamm = qammod(ham_encoded, 4, "gray", "InputType","bit","UnitAveragePower",true);
    sym_conv = qammod(conv_encoded, 4, "gray", "InputType","bit","UnitAveragePower",true);

    rec_hamm = awgn(sym_hamm, SNR(i));
    rec_conv = awgn(sym_conv, SNR(i));

    sym_hamm_rec = qamdemod(hamm_encoded, 4, "gray", "OutputType","bit");
    sym_conv_rec_hard = qamdemod(rec_conv, 4,"gray", "OutputType","bit");
    sym_conv_rec_soft = qamdemod(rec_conv, 4, "gray", "OutputType","llr");
end
