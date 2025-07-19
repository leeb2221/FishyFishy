clear, clc, close all

% Simulation type: 
% 1 - Particle Advection
% 2 - Velocity Heat Map
% 3 - Curl Heat Map
% 4 - Pressure Heat Map
TYPE = 3;

% Simulation parameters
s = 100;
ar = 5;
J = [0 1 0; 1 0 1; 0 1 0]/4;
n = 120000;

% Get the current jet colormap
JetMap = colormap("jet");
negJetMap = flipud(JetMap);

% Create a grid
[X, Y] = meshgrid(1:s*ar, 1:s);

% Initialize pressure and velocity fields
p = zeros(s, s*ar);
vx = zeros(s, s*ar);
vy = zeros(s, s*ar);

% Initial positions of particles
[px, py] = meshgrid(6:6, 1:s);

% Read the data from the Excel file
filename = 'fish.xlsx';
range = 'A1:DP100';
 
data1 = readmatrix(filename, 'Sheet', 1, 'Range', range);
mask1 = data1 ~= 1;
data2 = readmatrix(filename, 'Sheet', 2, 'Range', range);
mask2 = data2 ~= 1;
data3 = readmatrix(filename, 'Sheet', 3, 'Range', range);
mask3 = data3 ~= 1;
data4 = readmatrix(filename, 'Sheet', 4, 'Range', range);
mask4 = data4 ~= 1;
data5 = readmatrix(filename, 'Sheet', 5, 'Range', range);
mask5 = data5 ~= 1;
mask6 = flipud(mask2);
mask7 = flipud(mask3);
mask8 = flipud(mask4);
mask9 = flipud(mask5);

% The sequence of animation frames
sequence = [1 2 3 4 5 4 3 2 1 6 7 8 9 8 7 6];
numFramesPerMask = 3;  % change animation every n frames
FishCounter = 0;       % counts frames
masks = {mask1, mask2, mask3, mask4, mask5, mask6, mask7 mask8 mask9};


% Save these initial positions for the inflow
pxo = px;
pyo = py;
ramp = 9;
rampmax = 9;

mask = mask1;

f = figure(1);
set(f, 'WindowState', 'maximized');

% Main simulation loop (stops when closing figure window)
while ishandle(f)
    ramp = ramp+0.1;
    if ramp>rampmax
        ramp = rampmax;
    end
    start = tic;
    
    FishCounter = FishCounter + 1;

    % Determine which mask to use based on the frameCounter
    sequenceIndex = floor((FishCounter - 1) / numFramesPerMask);
    currentIndexInSequence = mod(sequenceIndex, length(sequence)) + 1;
    currentMaskIndex = sequence(currentIndexInSequence);
    mask = masks{currentMaskIndex};

    vx(1:1, :) = 0;
    vy(end:end, :) = 0;    

    % Compute right-hand side for pressure equation
    rhs = -divergence(vx, vy);

    % Jacobi iteration to solve for pressure
    for i = 1:1000
        % Update interior cells
        p = (conv2(p, J, 'same') + rhs/2);
        p(1:1,1:end) = p(2:2,1:end);                % Top Bound
        p(end:end,1:end) = p(end-1:end-1,1:end);    % Bottom Bound
        p(1:end,end:end) = 0;
        p(1:end,1:1) = ramp;
    end

    % Compute velocity gradient and update velocities for non-boundary pixels
    [dx, dy] = gradient(p);
    vx(2:end-1, 2:end-1) = vx(2:end-1, 2:end-1) - dx(2:end-1, 2:end-1);
    vy(2:end-1, 2:end-1) = vy(2:end-1, 2:end-1) - dy(2:end-1, 2:end-1);
    vx(~mask) = 0;
    vy(~mask) = 0;


    % Advect velocity field using Runge-Kutta 4th order method (-1 = backward)
    [pvx, pvy] = RK4(X, Y, vx, vy, -1);
    vx = interp2(vx, pvx, pvy, 'linear', 0);
    vy = interp2(vy, pvx, pvy, 'linear', 0);

    % Visualization of particle positions
    if TYPE == 1
        % Advect particles using Runge-Kutta 4th order method (1 = forward)
        [px, py] = RK4(px, py, vx, vy, 1);
    
        % Add the inflow particles
        px = [px; pxo];
        py = [py; pyo];
        % Remove particles after to many
        if length(px) > n
            px = px(end-n+1:end);
            py = py(end-n+1:end);
        end
        % Remove any outflow particles and ones that fall into circle
        index = (px < s*ar-2);
        px = px(index);
        py = py(index);
        
        % Plot mask as background
        imagesc(mask);
        colormap(gray);
        hold on;
    
        % Overlay particles
        scatter(px, py, 1, 'filled');
        hold off;
        axis equal;
        axis([0 s*ar 0 s]);
        xticks([]);
        yticks([]);
        
        
        stop = toc(start);
        FPS = 1/stop;
        title(sprintf("DIV/cell: %.5f    # Particles: %d    FPS: %.2f", ...
            sum(abs(divergence(vx(2:end-1,10:end-10), vy(2:end-1,10:end-10))), 'all')/(s*ar), length(px), FPS));
        drawnow;

    % Visualization of velocity field
    elseif TYPE == 2
        velocity_Mag = sqrt(vx(2:end-1,10:end).^2 + vy(2:end-1,10:end).^2);
        velocity_Mag = imresize(velocity_Mag, 10, 'bicubic');
        imagesc(velocity_Mag)
        colormap(JetMap)
        xticks([]);
        yticks([]);

        hold on;
    
        % Overlay current mask
        maskDisplay = double(~mask);
        maskDisplay = maskDisplay(2:end-1,10:end-10);
        maskDisplay = imresize(maskDisplay, 10, 'nearest'); 
        h = imagesc(maskDisplay);
        set(h, 'AlphaData', 1 * maskDisplay); % Transparency
        colormap(JetMap);
        hold off;
 
        stop = toc(start);
        FPS = 1/stop;
        title(sprintf("DIV/cell: %.5f FPS: %.2f", ...
            sum(abs(divergence(vx(2:end-1,10:end-10), vy(2:end-1,10:end-10))), 'all')/(s*ar), FPS));
        axis equal;
        drawnow;
    % Visualization of divergance field
    elseif TYPE==3
        DIVER = abs(curl(vx, vy));
        DIVER = DIVER(2:end-1,10:end-10);
        DIVER = imresize(DIVER, 10, 'bicubic');
        imagesc(DIVER)
        colormap(negJetMap)
        caxis([min(abs(curl(vx, vy)),[],"all") max(abs(curl(vx, vy)),[],"all")]);
%         axis([0 s*ar 0 s]);
        xticks([]);
        yticks([]);

        hold on;
    
        % Overlay current mask
        maskDisplay = double(~mask);
        maskDisplay = maskDisplay(2:end-1,10:end-10);
        maskDisplay = imresize(maskDisplay, 10, 'nearest'); 
        h = imagesc(maskDisplay);
        set(h, 'AlphaData', 1 * maskDisplay); % Transparency
        colormap(negJetMap);
        hold off;
    

        stop = toc(start);
        FPS = 1/stop;
        title(sprintf("DIV/cell: %.5f FPS: %.2f", ...
            sum(DIVER, 'all')/(s*ar), FPS));
        axis equal;
        drawnow;
    % Visualization of pressure field
    elseif TYPE == 4
        PRESS = p;
        PRESS = flipud(PRESS); % Flip for correct orientation
        PRESS = PRESS(2:end-1,10:end-10);
        PRESS = imresize(PRESS, 10, 'bicubic');
        imagesc(PRESS);
        colormap(JetMap);
        xticks([]);
        yticks([]);
        axis equal;
        drawnow;
    end
end

% Function for Runge-Kutta 4th order method for advection
function [x_new, y_new] = RK4(px, py, vx, vy, h)
   k1x = interp2(vx, px, py, 'linear', 0);
   k1y = interp2(vy, px, py, 'linear', 0);
   k2x = interp2(vx, px + h/2 * k1x, py + h/2 * k1y, 'linear', 0);
   k2y = interp2(vy, px + h/2 * k1x, py + h/2 * k1y, 'linear', 0);
   k3x = interp2(vx, px + h/2 * k2x, py + h/2 * k2y, 'linear', 0);
   k3y = interp2(vy, px + h/2 * k2x, py + h/2 * k2y, 'linear', 0);
   k4x = interp2(vx, px + h * k3x, py + h * k3y, 'linear', 0);
   k4y = interp2(vy, px + h * k3x, py + h * k3y, 'linear', 0);
   x_new = px + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
   y_new = py + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
end
