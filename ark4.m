%
%Adaptive fourth-order fixed-time-step Runge-Kutta ODE solver
%
%Parameters:
%   tZero   - Initial time (float)
%   deltaT  - Time step (float)
%   n       - Number of iterations (int)
%   xZero   - Initial x (vector)
%   fx      - ODE to solve (function)
%   tol     - Error tolerance (float)

function output = ark4(tZero, deltaT, n, xZero, fx, tol)
    
    output(:,1) = xZero;
    h = deltaT;
    t = tZero;
    
    for i = 1:n
        currPoint = output(:,i);
        
        while true
            tempOutput1 = rk(fx, h, t, currPoint);
            tempOutput2 = rk(fx, h/2, t, currPoint);
            tempOutput2 = rk(fx, h/2, t, tempOutput2);
            
            err = largestElem(tempOutput1, tempOutput2);
            
            if(err <= tol)
                break; 
                if(err <= tol/4)
                    h = h*2;
                end
            end
            
            h = h/2;
        end
    
        output(:,i+1) = tempOutput2;
        t = t+h;
    end
end

function tempOut = rk(fx, h, t, currPoint)

    k1 = h * fx(t, currPoint);
    k2 = h * fx(t+(t/2), currPoint+(k1/2));
    k3 = h * fx(t+(t/2), currPoint+(k2/2));
    k4 = h * fx(t, currPoint+k3);

    tempOut = currPoint + (1/6)*(k1+2*k2+2*k3+k4);
    
end

function err = largestElem(point1, point2)

    diff = abs(point1-point2);
    
    err = max(diff);
end