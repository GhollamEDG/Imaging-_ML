function ramp_break_idx  = sync_ramp(input_frame)
%% getting the signal ready
Rs = 2e6;        % sampling rate
N_spl = length(input_frame);
input_frame = input_frame - mean(input_frame);
% input_frame = input_frame / sqrt(sum(input_frame.^2)/length(input_frame));

% first diff
rs_dv = diff(input_frame.');

%% filter out the high freq "noise"
N   = 50;        
Fp  = 1e3;
Ap  = 0.01;
Ast = 80;

Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);
NUM = firceqrip(N,Fp/(Rs/2),[Rp Rst],'passedge');
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
        [mininum, ind] = min(input_frame(indices(i)-threshold:indices(i)));
    else
        [maximum, ind] = max(input_frame(indices(i)-threshold:indices(i)));
    end
    indices(i) = ind + indices(i) - threshold -1;
end
IDX = zeros(1,length(input_frame));
IDX(indices) = input_frame(indices);

ramp_break_idx = find(IDX ~=0);

% plot the results
              
% A = 0.6e2; % normalizing factor
% figure;
% plot(1:N_spl, input_frame, '.', 1:N_spl, A*[rs_dv,0],  ...
%     1:N_spl, A*[rs_dv_lp,0], '.', 1:N_spl, A^2*rs_dv2, '*', ...
%     1:N_spl, IDX, '*');
% legend('original signal', 'differenced once',...
%     'lowpass + diff1', 'diff + lowpass + diff', 'predicted');
% 
figure;
plot(1:N_spl, input_frame, '.', ...
    1:N_spl, IDX, '*');
legend('original signal', 'predicted');


end