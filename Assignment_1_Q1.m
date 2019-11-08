
convolute([1,3,2,1,2,2,5],[1,0,1,0,1])

function conv = convolute(x1, x2)

    % assign smaller array to be h, then flip it.
    if length(x1) > length(x2)
        x = x1;
        h = flip(x2);
    else
        x = x2;
        h = flip(x1);
    end
    % stores the original size of x for iteration
    x_size = length(x);
    % pad x with zeroes of length h
    x = [zeros(1, length(h) - 1) x zeros(1, length(h) - 1)];
    % instantiate empty output vector
    y = int16.empty;
    
    % for each element in x, multiply it by every element in h and append
    % the sum to y.
    for i = 1:x_size + length(h) - 1
        arg1 = x(i);
        sum = 0;
           for j = 1:length(h)
               sum = sum + h(j) * x(j + i - 1)  ; 
           end
        y(i) = sum;
    end
    conv = y;
end