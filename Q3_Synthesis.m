% Assignment 1 Question 3
% Systems Design Engineering


function [Xt] = MyiFT( Xw, w, t )
      
    e  = @(t,w) exp(-1j*t*w);
    
    Xt = zeros(1,length(t)) % create an array of zeros that will eventually store values of the inverse fourier transform function
    
    for i = 1:length(t) %calculating result at each frequency
        result = 0; % store result of integral at given Xw
        
        for j = 1:length(w) % looping through each value of t
            result = result + ( Xw(j) * e(w(j), t(i)) );
        end
             
        Xt(i) = (1/2*pi) * result; % stores the result of fourier transform at each respective t
    end
end

