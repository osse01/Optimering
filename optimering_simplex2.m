clear all
clc
close all

maxIter = 20;

k = 5;
f1 = @(x1,x2) 10^k*(x2-x1.^2).^2 + (1-x1).^2;
x = [[-5,5]',[5,2]',[-5,-5]']

figure
axis(-5,5,-5,5)
hold on
gifFile = 'simplex.gif';
gradient = logspace(0,log10(1/maxIter), maxIter) % (logarithmic) Gradient scale (0 to 1)
for i = 1:maxIter
    x = simplex(f1,x);
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
