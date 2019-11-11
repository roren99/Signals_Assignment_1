 
sample_rate = 0.01
t = 0.0:sample_rate:3;
xt = exp(-j*t); %sin(8*t); %1+cos(2*pi*t)/4 + cos(2*pi*t*2)/2 + cos(2*pi*t*3)/3;

% 1+cos(2*pi*t)/4 + cos(2*pi*t*2)/2 + cos(2*pi*t*3)/3 in frequency domain:
ak = zeros(1,41);
ak(21) = 1;
ak(20)= 0.1134;
ak(22)= 0.1134;
ak(19)= 0.2389;
ak(23)= 0.2389;
ak(18)= 0.163;
ak(24)= 0.163;

[ak, k] = fourier_series(xt, t);
[retrieved_xt, retreived_t] = inv_fourier_series(ak, k);

function [dv, iv] = fourier_series(xt, t)
   figure('Name', 'Inputted x(t)'), plot(t,xt);
   dt = t(2) - t(1);
   [period, N] = get_ft(xt, t, log10(1/dt) - 1);
   omega = (2 * pi) / period;
   
   integration_interval = dt : dt : period;
   k_min = -1 * 20;
   k_max = 20;
   ks = k_min:1:k_max;
   ak = zeros(1, length(ks)); % instantiate array of aks as zeroes

   xt_range = xt(1:N)
   
   figure('Name', 'Fundamental Period'), plot(integration_interval, xt_range)

   s = sum((xt(1: N) .* exp(-1i * omega * 1 * integration_interval ) * dt))
   
   for k = ks
       integral = sum((xt(1: N) .* exp(-1i * omega * k * integration_interval ) * dt));
       ak(k + k_max + 1) = integral / period;
   end
   dv = ak;
   iv = ks;
   figure('Name', 'ak real'), plot(ks,real(ak));
   figure('Name', 'ak imaginary'), plot(ks,imag(ak));
end

function [dv, iv] = inv_fourier_series(ak, ks)
    ts = 0:0.01:3;
    omega = 2 * pi;
    xts = []; % instantiate array of xts
    for singular_t = ts
        xts_sum_for_one_t = sum(ak .* exp(-1i * omega * singular_t * ks));
        xts = [xts xts_sum_for_one_t];  
    end
    dv = xts;
    iv = ts;
    figure('Name', 'retreived X(t) vs t'), plot(ts,xts);
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
        
        if (all(abs(left_period - right_period) < TOLERANCE) == ones(1, length(left_period)))
            period = t(length(left_period) + 1) - t(1);
            period_index_count = length(left_period);
        end
    end

end