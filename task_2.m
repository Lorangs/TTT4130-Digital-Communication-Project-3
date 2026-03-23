% Convolutional code

k = 3;
msg_len = 10;
ramp_down = true;		% will ramp down the convolutional encoder at the end by padding the input message with trailing k-1 zeros.
% generate test message
msg = randi([0 1], 1, msg_len);
disp("Initial message");
disp(msg);

% Generator vector for the shifted values of msg
g1 = [1 1 1];
g2 = [1 0 1];

function encoded_msg = encode(msg, k, g1, g2, ramp_down) 
	if ramp_down
		padded_msg = [zeros(1, k-1) msg zeros(1, k-1)];
	else 
		padded_msg = [zeros(1 ,k-1) msg];
	end

	shift_reg = zeros(1, k);
	encoded_msg = zeros(1, 2*length(msg) + 2*(k-1));

	for i = 1:length(padded_msg)-k+1
		shift_reg = padded_msg(i:i+k-1);				% update the shift register
		encoded_msg(2*i-1) 	= mod(shift_reg * g1.', 2);	% output of generator branch 1
		encoded_msg(2*i) 	= mod(shift_reg * g2.', 2);	% output of generator branch 2
	end
end 

encoded = encode(msg, k, g1, g2);
disp("encoded message");
disp(encoded);


% compare with built in matlab function
trellis = poly2trellis(k, [7 5]);

convencoder = comm.ConvolutionalEncoder(trellis);

if ramp_down
	encoded_msg_Comm = convencoder([msg, zeros(1,k-1)].');
else 
	encoded_msg_Comm = convencoder(msg.');
end
disp("test with matlab encoder");
disp(encoded_msg_Comm.');

exit;

