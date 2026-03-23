% Implement the Convolutional Decoder and encode a binary message 

m = 4; 

msg = [1 0 1 0]; % m binary digits
K = 3;

% Connection vectors 
g1 = [1 1 1];
g2 = [1 0 1];

pad = zeros(1,K-1);
padded_msg = [pad msg];

encoded = [];


for i = 1:length(padded_msg) - K + 1 
    v1 = mod(sum(padded_msg(i:i+K-1) .* g1), 2);
    v2 = mod(sum(padded_msg(i:i+K-1) .* g2), 2);
    encoded = [encoded, v1, v2];
end




% Encode using the built-in convolutional encoder for comparison
trellis = poly2trellis(3, [7 5]); % 7 is the binary code for g1, and 5 is binary code for g2


convencoder = comm.ConvolutionalEncoder(trellis);
codeword = convencoder(msg.');

% Display for comparison
disp('Custom Encoded Message:');
disp(encoded);
disp('Built-in Encoder Codeword:');
disp(codeword.');