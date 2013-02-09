function rk4
    output = ruk4(0, 0.005, 10000, [3; 0.1], @dampedPendulum);
    figure('Name', 'x1 and x2');
    plot(output(1,:));
    hold on
    plot(output(2,:), 'color', 'green')
    hold off
end

function output = ruk4(tZero, deltaT, n, xZero, fx)
    
    output(:,1) = xZero;
    h = deltaT;
    %t = tZero:h:n/h;
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

%Function that calculates the motion of a forced, damped pendulum
function xprime = dampedPendulum(t, x)
%Static values
m = 0.1;
l = 0.1;
beta = 0;
alpha = 0;
A = 0;
g = 9.81;

xprime = [x(2); (A*cos(alpha*t) - beta*l*x(2) - m*g*sin(x(1)))/m*l];
end