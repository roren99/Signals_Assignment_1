 
sample_rate = 0.01
t = 0.1:sample_rate:3;
xt = 1+cos(2*pi*t)/4 + cos(2*pi*t*2)/2 + cos(2*pi*t*3)/3;

[ak, k] = fourier_transform(xt, t);
% [xt_2, t_2] = fourier_transform(xt, t);

function [dv, iv] = fourier_transform(xt, t)
   fig = figure, plot(t,xt);
   dt = t(2) - t(1);
   [period, period_index_count] = get_ft(xt, t, log10(1/dt) - 1);
   omega = (2*pi)/period;
   
   integration_range = 0 : dt : period - dt;
   
   k_min = -20;
   k_max = 20;
   ks = k_min:1:k_max;
   ak = zeros(1, length(ks)); % instantiate array of aks as zeroes

   for k = ks
       integral = sum((xt(1: period_index_count) .* exp(-1i * omega * k * integration_range ) * dt));
       ak(k + 20 + 1) = integral / period;
   end
   dv = ak;
   iv = ks;
   figure, plot(ks,real(ak));
end

function [dv, iv] = inv_fourier_transform(ak, k)

    ts = 0.1:0.1:3;
    xts = []; % instantiate array of xts
    
    for singular_t = ts
        xts_sum_for_one_t = sum(ak .* exp(-1i * omega * singular_t * k));
        xts = [xts xts_sum_for_one_t];
    end
        

end

function [period, period_index_count] = get_ft(xt, t, decimal_palces)
    MIN_PERIOD_SAMPLES = 3;
    TOLERANCE = 0.15;
    sum_of_first_two_outputs = xt(1) + xt(2);
    decimal_places = 2
    sum = xt(3);
    
    left_period = xt(1 : length(xt)/2);
    right_period = xt(length(xt)/2 : length(xt)/2 + length(left_period) - 1);
    
    period = -1; %default values
    period_index_count = 0;
    
    while length(left_period) > MIN_PERIOD_SAMPLES
        right_element_of_left_period = left_period(end);
        left_period(end) = [];
        right_period = [right_element_of_left_period right_period];
        right_period = right_period(1: length(right_period) - 2);
        
        a = left_period - right_period;
        b = abs(left_period - right_period) < TOLERANCE;
        c = all(abs(left_period - right_period) < TOLERANCE);
        
        if (all(abs(left_period - right_period) < TOLERANCE) == ones(1, length(left_period)))
            period = t(length(left_period) + 1) - t(1);
            period_index_count = length(left_period);
        end
    end

end