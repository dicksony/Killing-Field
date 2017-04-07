%%Tests for the BH function

R = 1; %Particle Initial radius
E = 1; %Particle Initial Energy
L = 1; %Particle angular momentum
Type = 0; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 30000; %Length of returned coordinate arrays
time_step = 0.01;
MASSIVE = 1;
ingoing = 1;

[radius, angle, time] = BH(Type,BH_Ang,BH_Mass,R,E,L,MASSIVE,SIZE,time_step,ingoing);

if ~isempty(time)
    polarplot(angle,radius); %Display data
end

%%Good Newtonian orbit Demo, elliptic
%{
R = 10; %Particle Initial radius
E = 0.55; %Particle Initial Energy
L = sqrt(120); %Particle angular momentum
Type = 0; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 30000; %Length of returned coordinate arrays
time_step = 0.01;
MASSIVE = 1;
ingoing = 1;
%}

%%Newtonian Circular orbit
%{
R = 1; %Particle Initial radius
E = 1; %Particle Initial Energy
L = 1; %Particle angular momentum
Type = 0; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 30000; %Length of returned coordinate arrays
time_step = 0.01;
MASSIVE = 1;
ingoing = 1;
%}