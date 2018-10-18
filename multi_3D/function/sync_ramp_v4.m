function idx = sync_ramp_v4(input_frame, ifplot)

N_spl = length(input_frame);
input_frame(input_frame > 0) = 0.5* input_frame(input_frame > 0);
input_frame = transpose(input_frame) - mean(input_frame);
input_frame = input_frame / vecnorm(input_frame);
input_frame_abs = abs(input_frame);
maxabs = max(input_frame_abs);
threshold = 0.97;
I = find(input_frame_abs > threshold *maxabs);
indices = [I(1), I(find(diff(I) > 5)+1)]+1;
I2 = [];
for i = 1: 2: length(indices)-1
    potential_point = floor(mean(indices(i:i+1)));
    if(input_frame(potential_point) > input_frame(indices(i)))
        [~, new_idx] = max(input_frame(indices(i):indices(i+1)));
    else
        [~, new_idx] = min(input_frame(indices(i):indices(i+1)));
    end
    new_idx = new_idx + indices(i) - 1;
    I2 = [I2, new_idx];
end

if(isempty(I2))
    idx = 0;
end
if(length(I2) == 1)
    idx = I2(1);
end
if(length(I2) >=2)
    if(I2(2) - I2(1)) > length(input_frame)/2
        idx = I2(2);
    else
        idx = I2(1);
    end
end



if ifplot
    
    IDX = zeros(1,N_spl);
    IDX(I2) = input_frame(I2);

    figure
    plot(1:N_spl, input_frame, 1:N_spl, IDX, '*');
    legend('original signal' , 'predicted');
end




end