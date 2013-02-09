function main
    output = ruk4(0, 0.005, 10000, [3; 0.1], @dampedPendulum);
    
    figure('Name', 'x1 and x2');
    plot(output(1,:));
    hold on
    plot(output(2,:), 'color', 'green')
    hold off
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