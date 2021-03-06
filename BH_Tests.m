%%Tests for the BH function

%%JAYMY UZE MII
R = 9.606; %Particle Initial radius
E = 0.5; %Particle Initial Energy
L = 4.2; %Particle angular momentum
Type = 0; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = -0.1; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 1000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 0;
ingoing = 1;

[radius, angle, time] = BH(Type,BH_Ang,BH_Mass,R,E,L,MASSIVE,SIZE,time_step,ingoing);

if ~isempty(time)
    polarplot(angle,radius); %Display data
    hold on
end

R = 9.606; %Particle Initial radius
E = 0.5; %Particle Initial Energy
L = 4.2; %Particle angular momentum
Type = 1; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = -0.1; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 1000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 0;
ingoing = 1;

[radius, angle, time] = BH(Type,BH_Ang,BH_Mass,R,E,L,MASSIVE,SIZE,time_step,ingoing);

if ~isempty(time)
    polarplot(angle,radius); %Display data
    hold on
end

R = 9.606; %Particle Initial radius
E = 0.5; %Particle Initial Energy
L = 5; %Particle angular momentum
Type = 2; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = 0.01; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 1000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 0;
ingoing = 1;

[radius, angle, time] = BH(Type,BH_Ang,BH_Mass,R,E,L,MASSIVE,SIZE,time_step,ingoing);

if ~isempty(time)
    polarplot(angle,radius); %Display data
    hold on
end

%%STANDARD TESTS
%{
R = 9.606; %Particle Initial radius
E = 1.001; %Particle Initial Energy
L = 4; %Particle angular momentum
Type = 2; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = 0.05; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 7000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 1;
ingoing = 1;

[radius, angle, time] = BH(Type,BH_Ang,BH_Mass,R,E,L,MASSIVE,SIZE,time_step,ingoing);

if ~isempty(time)
    polarplot(angle,radius); %Display data
end

%}

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

%%Hilarious Newtonian Ellipse
%{
R = 10; %Particle Initial radius
E = 0; %Particle Initial Energy
L = 1.4; %Particle angular momentum
Type = 0; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 30000; %Length of returned coordinate arrays
time_step = 0.03;
MASSIVE = 1;
ingoing = 1;
%}

%%Schwarzchild massive particle deflection
%{
R = 50; %Particle Initial radius
E = 7; %Particle Initial Energy
L = 80; %Particle angular momentum
Type = 1; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 700; %Length of returned coordinate arrays
time_step = 0.03;
MASSIVE = 1;
ingoing = 1;
%}

%%Schwarzchild MAssive particle time contraction
%{
R = 40; %Particle Initial radius
E = 1; %Particle Initial Energy
L = 4; %Particle angular momentum
Type = 1; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 25000; %Length of returned coordinate arrays
time_step = 0.5;
MASSIVE = 1;
ingoing = 1;
%}

%%Schwarzschild Massive particle bowtie
%{
R = 40; %Particle Initial radius
E = 1; %Particle Initial Energy
L = 4.2; %Particle angular momentum
Type = 1; % 0-Newton/1-Scwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 20000; %Length of returned coordinate arrays
time_step = 0.03;
MASSIVE = 1;
ingoing = 1;
%}

%%Schwarzchild TRUE BEAUTY DEFINED
%{
R = 9.606; %Particle Initial radius
E = sqrt(2*(-0.0214)+1); %Particle Initial Energy
L = 4; %Particle angular momentum
Type = 1; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = 0; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 50000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 1;
ingoing = 1;
%}

%%Kerr Counter-rotation Demo
%{
R = 9.606; %Particle Initial radius
E = 1.001; %Particle Initial Energy
L = 4; %Particle angular momentum
Type = 2; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = 0.05; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 7000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 1;
ingoing = 1;
%}

%%Kerr Direct comparison with above, co-rotation
%{
R = 9.606; %Particle Initial radius
E = 1.001; %Particle Initial Energy
L = 4; %Particle angular momentum
Type = 2; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = -0.05; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 7000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 1;
ingoing = 1;
%}

%%Schwarz Light BENDS?!?!?!?!
%{
R = 9.606; %Particle Initial radius
E = 0.5; %Particle Initial Energy
L = 4.2; %Particle angular momentum
Type = 1; % 0-Newton/1-Schwarz/2-Kerr
BH_Ang = -0.1; %Black Hole Angular Momentum
BH_Mass = 1; %Black Hole Mass
SIZE = 1000; %Length of returned coordinate arrays
time_step = 0.08;
MASSIVE = 0;
ingoing = 1;
%}