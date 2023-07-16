%Project 2 Hitting a home run, with air resistance
%Name:Arvind Purohit
%Description:Exporting data and analyzing it in Excel
clear
clf 

% ----- define given information -----

x0 = 0; y0 = 0; %Position

v0mph = 112;    % velocity, in mph 
phi0deg = 32;   % launch angle, in deg
g = 10;         % gravity, N/kg 
m = 0.142;      % mass of a baseball, in kg
mtoft = 3.28084; %Conversion factor
% ----- set up more variables -----

mph2mps = 5280*12*2.54/100/3600; %mph to m/s conversion
deg2rad = pi()/180;   % degrees to radians conversion
 
v0 = v0mph*mph2mps;    % exit velocity, in m/s 
phi0 = phi0deg*deg2rad;  % launch angle, in rad 
 
v0x = v0*cos(phi0);  % x-component of initial velocity in m/s
v0y = v0*sin(phi0);  % y-component of initial velocity in m/s

t_H =  v0y/g; % time to reach max. height H, in s
T = 2*t_H; % time to land on ground, in s
H = t_H*v0y/2; % maximum height, in m
R = v0x*T; % range, in m
 
% ----- set up a time array and compute x(t) and y(t) analytically -----

tmin = 0; 
tmax = T; % based on the t eq written above 
N = 2000; % intervals 
t = linspace(tmin, tmax, 1+N); % time array, connects x(t) and y(t)
 
xt = x0 + v0x*t; 
xf = xt*mtoft; % converts xt to ft 
yt = y0 + v0y*t - (1/2)*g*t.^2; 
yf = yt*mtoft; % convert yt to ft 

dt = (tmax-tmin)/N; % step size
y = zeros(1, 1+N);  % initialize y(t)
x = zeros(1, 1+N);  % initialize x(t)

y(1) = y0; 
x(1) = x0; 

vy = v0y; 
vx = v0x; 

%---Components for the Drag -----
C = input("Choose a value for C (Between 0.2-0.3): " );
A = 0.004301; %in m^2
Pa = 1.225; %in Kg/m^3
con = 0.5*C*Pa*A;
for n = 1:N 
    v = sqrt(vx^2 +vy^2); %Speed of baseball, in m/s
    Fy = -m*g - con*v*vy; % force in y direction due to gravity and drag.
    ay = Fy/m; % acceleration in the y direction
     
    Fx = - con*v*vx; % force in the x-acceleration due to drag 
    ax = Fx/m; % acceleration in the x-direction 
    
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;
    vy = vy + ay*dt; 
    
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2; 
    vx = vx + ax*dt;
    
end
 
y = y*mtoft; %converts y to ft
x = x*mtoft; %converts x to ft 

%Plotting

p = plot(xf, yf, x, y,'LineWidth',4);
grid on
grid minor

ax = gca; 
ax.FontSize = 16; 
ax.GridAlpha = 0.4;
ax.MinorGridAlpha = 0.5;

xlabel('x (ft)', 'FontSize', 18) 
ylabel('y (ft)', 'FontSize', 18)

title({'ECE 202, Project 2, Phase 3: Trajectory of a baseball', ...
    'With drag Vs. No Drag'},'FontSize', 24)
l1 = sprintf('With Drag, C = %g',C);
legend('Without Drag',l1,'FontSize', 18)

ylim([-5,120])

%Check
 
CheckSumy = sum(abs(y-yf)) % Check should be zero 
CheckSumx = sum(abs(x-xf)) % Check should be zero 

%----all elements

mat = [t; x; y;].';
writematrix(mat, "val.csv")