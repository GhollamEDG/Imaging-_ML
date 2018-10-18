function high_freq_start = sync_ramp(input_frame,Rs)

% sampling rate\
N_spl = length(input_frame);
input_frame = transpose(input_frame) - mean(input_frame);
smaples = input_frame / vecnorm(input_frame);
level = 1;
[c,l] = wavedec(input_frame,level,'haar');

d1 = detcoef(c,l,level);
rs_dv = -interpft(d1,2*length(d1));

% %% filter out the high freq "noise"
% N   = 50;
% Fp  = 1e3;
% Ap  = 0.01;
% Ast = 80;
% 
% Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
% Rst = 10^(-Ast/20);
% NUM = firceqrip(N,Fp/(Rs/2),[Rp Rst],'passedge');
% rs_dv_lp = filter(NUM, 1, rs_dv);

%% second differentiation after filtering
level = 1;
[c,l] = wavedec(rs_dv,level,'haar');

d1 = detcoef(c,l,level);
rs_dv2 = -interpft(d1,2*length(d1));

%% find the peaks

I = find(abs(rs_dv2) > 0.9*max(abs(rs_dv2)));
indices = [I(1), I(find(diff(I) > 100)+1)]+1;
while(indices(1) < 10)
    indices = indices(2:end);
end
% figure
% plot(1:N_spl, input_frame, 1:N_spl, 10*rs_dv);
% figure;
% plot(1:length(abs(rs_dv2)), abs(rs_dv2));
for i=1:1:length(indices) % revise the indices
    threshold = min(indices(i)-1, 5);
    if (rs_dv2(indices(i)) < 0)
        [~, ind] = max(input_frame(indices(i):indices(i)+threshold));
    else
        [~, ind] = min(input_frame(indices(i):indices(i)+threshold));
    end
    indices(i) = ind + indices(i) -1;
end
 %% decide what to return
 


IDX = zeros(1,N_spl);
IDX(indices) = input_frame(indices);


if(isempty(indices))
    high_freq_start = 0;
end
if(length(indices) == 1)
    high_freq_start = indices(1);
end
if(length(indices) >=2)
    if(indices(2) - indices(1)) > (Rs*6e-3/2)
        high_freq_start = indices(2);
    else
        high_freq_start = indices(1);
    end
end
% A = 0.6e1;
% figure
% plot(1:N_spl, input_frame, '.', 1:N_spl, A*rs_dv, 1:N_spl, IDX, '*');
% legend('original signal', 'differenced once', 'predicted');
end
