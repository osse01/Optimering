clear all
clc
close all

maxIter = 25;

k = 1;
f1 = @(x1,x2) 10^k*(x2-x1.^2).^2 + (1-x1).^2;
%x = [[0,5]', [0,0]', [5,5]'] % 1
%x = [[2,-3]', [4,2]',[-5,5]'] % 2

f2 = @(x1, x2) -4.*exp(-( (x1+2).^2 + (x2+1).^2 )/10 ) + ...
                4.*exp(-( (x1+2).^2 + (x2+1).^2 )/100) + ...
                1 / 100*( (x1+2).^2 + (x2+1).^2 + x1);
%x = [[10,10]', [-7,5]', [0,-6]'] % 3
%x = [[5,10]', [-7,5]', [0,-10]'] % 4

figure
axis([-5,5,-5,5]) % For f1
%axis([-10,10,-10,10]) % For f2
f = f1;

hold on
gifFile = 'simplex.gif';
gradient = logspace(0,log10(1/maxIter), maxIter); % (logarithmic) Gradient scale (0 to 1)

plot( [x(1,1:end) x(1,1)],[x(2,1:end) x(2,1)], ...
          'Color', 'black', 'LineWidth', 2);
f_middle = zeros(1,maxIter);
for i = 1:maxIter
    x = simplex(f,x);
    
    % Rest is for plotting
    middle = [1/3*sum(x(1,1:end))  1/3*sum(x(2,1:end))]
    f_middle(i) = f(middle(1), middle(2));
    color_segment = [1-gradient(i),0,gradient(i)];
    plot( [x(1,1:end) x(1,1)],[x(2,1:end) x(2,1)], ...
          'Color', color_segment, 'LineWidth', 2);


    % Capture the frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if i == 1
        % First frame: create the file
        imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        % Subsequent frames: append to the file
        imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
% Plot the middle point function value
figure
plot(f_middle);
