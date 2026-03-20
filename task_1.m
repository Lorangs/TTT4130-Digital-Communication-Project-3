% Hamming code

k = 4;
n = 7;
m = n - k;

num_errors = 3;

% Generate different seeds for random functions
rng('shuffle');

% example of known message
x = randi([0 1], 1, k);
disp("x:"); 
disp(x);

[H, G] = hammgen(m);
code = encode(x, n, k);
disp("code");
disp(code);

% add parity as MSB of code
parity_bit = mod(sum(code), 2);
code = [parity_bit, code];
disp("Parity-bit");
disp(parity_bit);
disp("code with parity (MSB)");
disp(code);

% Syndrome table
LUT = syndtable(H);
%disp("LUT:");
%disp(LUT);

SyndTable = mod(LUT * H.', 2);
%disp("Syndrome Table:");
%disp(SyndTable);

% generate a bit flip at random location
e = zeros(1, n+1); % plus paritybit
e_idx = randperm(7, num_errors);
e(e_idx) = 1;
disp("Error vector");
disp(e);


r = mod(code + e, 2);
disp("Received Signal");
disp(r);

% get received parity bit
r_p = r(1);

% calculate paritybit of received signal
est_p = mod(sum(r(2:end)), 2);

% check if parity bit is flipped. 
parity_check_n = xor(est_p, r_p);
disp("ParityCheck_N");
disp(parity_check_n);

% calculate syndrom of received signal
s = mod(r(2:end) * H.', 2);
disp("Syndrome:");
disp(s);

[~, calc_e_idx] = ismember(s, SyndTable, 'rows');
disp("Error Index at");
disp(calc_e_idx);


recoverable = true;
if parity_check_n == 0
	if all(s==0)
		disp("No errors detected");
	else
		disp("Two bit errors detected. Not recoverable");
		recoverable = false;
	end
else
	if all(s==0)
		disp("Parity bit flipped.");
		r(1) = mod(r(1) + 1, 2);
	else
		disp("Single bit error detected and corrected. There is no guarantee there are 3 or more errors lurking.");
		r(2:end) = mod(r(2:end) + LUT(calc_e_idx, :), 2);
	end
end

if recoverable
	decoded_msg = decode(r(2:end), n, k);
	disp("Decoded bits");
	disp(decoded_msg);
end
exit;
