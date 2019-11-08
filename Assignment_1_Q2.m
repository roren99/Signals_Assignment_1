% t = 0.1:1:3;
% xt = 1+cos(2*pi*t)/4 + cos(2*pi*t*2)/2 + cos(2*pi*t*3)/3;
t = 0:0.1:3;
fourier_transform(cos(2*pi*t), t);

function ctft = fourier_transform(xt, t)
   fig = figure, plot(t,xt);
   period = get_ft(xt, t);
   omega = (2*pi)/period;
   
   dt = 0.01;
   integration_range = (period * -1) / 2 : dt : period / 2;
   
   ks = 0:1:(length(integration_range) - 1);
   ak = zeros(1, length(ks)) % instantiate array of aks as zeroes
   
   for k = 1:length(integration_range)      
       integral = 0;
       for b = 1:length(integration_range)
           integral = integral + (xt .* exp(-1i * omega * k * t ) * dt);
       end
       ak(k) = sum(integral) / period;
   end
    figure, plot(ks,ak);
end

function period = get_ft(xt, t)
    TOLERANCE = 0.01;
    sum_of_first_two_outputs = xt(1) + xt(2);
    
    sum = xt(3);
    
    %the period occurs when the sum of two consecutive elements are within
    %some tolerance of eachother.
    for i = 4:length(xt)
        sum = sum + xt(i);
        if (sum < sum_of_first_two_outputs + TOLERANCE) && (sum > sum_of_first_two_outputs - TOLERANCE)
            %the actual period is the delta x values corresponding with the
            %matching xt values.
            period = t(i) - t(1);
            return
        end
        sum = sum - xt(i - 1);
    end
    period = -1;
end