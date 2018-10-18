%% load the signal and initiations
close all; clear; clc;
samples = transpose(read_complex_binary2('../data/0918_m_d1_4_2.dat')); % read the samples from file
% since reading the samples from file take considerable time, you may
% comment this line after running it once.

period = 12e3;              % the signal is periodic with a period of 12000 samples

fs = 2e6;        % sampling rate

%% getting the signal ready
total_length = length(samples);
rows = floor(total_length/period);
actual_length = period * rows;
samples = samples(1:actual_length);
samples = real(samples);
samples = samples - mean(samples);
% smaples = samples / sqrt(sum(samples.^2)/length(samples));
% first diff
rs = samples(period+1+fs:2*period+fs); 
rs_dv = diff(rs);

%% filter out the high freq "noise"
N   = 50;        
Fp  = 1e3;
Ap  = 0.01;
Ast = 80;

Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);
NUM = firceqrip(N,Fp/(fs/2),[Rp Rst],'passedge');
rs_dv_lp = filter(NUM, 1, rs_dv);

%% second differentiation after filtering
level = 1;
[c,l] = wavedec(rs_dv_lp,level,'haar');

d1 = detcoef(c,l,level);
rs_dv2 = -interpft(d1,2*length(d1));

%% find the peaks

I = find(abs(rs_dv2) > 1.4e-5);
threshold = 30;
indices = [I(1), I(find(diff(I) > 500)+2)];
for i=1:1:length(indices) % revise the indices
    if(rs_dv2(indices(i)) > 0)
        [mininum, ind] = min(rs(indices(i)-threshold:indices(i)));
    else
        [maximum, ind] = max(rs(indices(i)-threshold:indices(i)));
    end
    indices(i) = ind + indices(i) - threshold -1;
end
IDX = zeros(1,period);
IDX(indices) = rs(indices);

%% plot the results
              
A = 0.6e2; % normalizing factor
figure;
% plot(1:period, 3e2*[rs_dv,0], 1:period, rs, '.');
plot(1:period, rs, '.', 1:period, A*[rs_dv,0],  ...
    1:period, A*[rs_dv_lp,0], '.', 1:period, A^2*rs_dv2, '*', ...
    1:period, IDX, '*');
legend('original signal', 'differenced once',...
    'lowpass + diff1', 'diff + lowpass + diff', 'predicted');

figure;
% plot(1:period, 3e2*[rs_dv,0], 1:period, rs, '.');
plot(1:period, rs, '.', ...
    1:period, IDX, '*');
legend('original signal', 'predicted');

% figure;
% plot(1:period, [rs_dv,0]);