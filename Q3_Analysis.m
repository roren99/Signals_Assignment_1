% Assignment 1 Question 3
% Systems Design Engineering

t = 0.1:0.1:8;
Xt = cos(t)
w = 0:0.1:4

function [Xw] = MyFT(Xt, t, w)

    e  = @(t,w) exp(-1j*t*w);
    
    Xw = zeros(1,length(w)) % create an array of zeros that will eventually store values of the fourier transform function
    
    for i = 1:length(w) %calculating result at each frequency
        result = 0; % store result of integral at given Xt
        
        for j = 1:length(t) % looping through each value of t
            result = result + ( Xt(j) * e(t(j), w(i)) );
        end
             
        Xw(i) = result; % stores the result of fourier transform at each respective w
    end
    plot(w, Xw);
end
