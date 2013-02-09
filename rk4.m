%
%Fourth-order fixed-time-step Runge-Kutta ODE solver
%
%Parameters:
%   tZero   - Initial time (float)
%   deltaT  - Time step (float)
%   n       - Number of iterations (int)
%   xZero   - Initial x (vector)
%   fx      - ODE to solve (function)

function output = rk4(tZero, deltaT, n, xZero, fx)
    
    output(:,1) = xZero;
    h = deltaT;
    t = tZero;
    
    for i = 1:n
        
        k1 = h * fx(t, output(:,i));
        k2 = h * fx(t+(t/2), output(:,i)+(k1/2));
        k3 = h * fx(t+(t/2), output(:,i)+(k2/2));
        k4 = h * fx(t, output(:,i)+k3);
    
        output(:,i+1) = output(:,i) + (1/6)*(k1+2*k2+2*k3+k4);
        t = t+h;
    end
end
