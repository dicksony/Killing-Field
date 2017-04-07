% BH
%   Args:
%   BH_Type: 0/1/2 = Newton/Schwarz/Kerr
%   BH_L: BH Angular momentum/mass, 0<=BH_L<1, only used in Kerr case
%   BH_M: BH Mass
%   P_R: Particle initial radius
%   P_E: Particle initial energy, dependent on BH Type:
%     Newtonian: Radial kinetic energy
%     Schwarz/Kerr: Energy per unit mass at infinity, (P_E^2-1) ~ Kinetic_E
%   P_L: Particle specific (per unit mass) angular momentum 
%   P_MASSIVE: 1 for massive particle, 0 for photon
%   POSITION_ARRAY_SIZE: REquested length of coordinate arrays
%   d_tau: Requested time step size in proper time (Reflects accuracy)
%   ingoing_flag: 1 if particle starts infalling or 0 for outgoing
%
% NOTE THAT GEOMETRIZED UNITS ARE USED THROUGHOUT (G=1,c=1)
% Returns empty array for pos_t if invalid parameters chosen

%%TODO: FIX THE NEWTONIAN HYPERBOLIC/

function [pos_r,pos_phi,pos_t] = ...
    BH(BH_Type,BH_L,BH_M,P_R,P_E,P_L,P_MASSIVE,POSITION_ARRAY_SIZE,d_tau,...
    ingoing_flag)

%%TODO: Add support for massless particles.

pos_r = zeros(1,POSITION_ARRAY_SIZE);
pos_r(1) = P_R;
pos_t = zeros(1,POSITION_ARRAY_SIZE);
pos_phi = zeros(1,POSITION_ARRAY_SIZE);

if P_MASSIVE == 1
    K = 1;
else
    K = 0;
end

if BH_Type == 1
    
    %Check for valid energy/Angular momentum values
    if (P_E^2 - (1-2*BH_M/pos_r(1))*(1+P_L^2/pos_r(1)^2)) < 0
        fprintf('Invalid Energy/Angular Momentum Selection');
        pos_t = [];
        return;
    end
    
    %Find Turning points based on tailored potential
    potentialPolynom = [(P_E^2-1) 2*BH_M -P_L^2 2*BH_M*P_L^2];
    turningRoots = roots(potentialPolynom);
    turningPoints = zeros(1,3);
    for i = 1:length(turningRoots)
        if turningRoots(i) > 0 && imag(turningRoots(i)) == 0
            turningPoints(i) = turningRoots(i);
        end
    end
    turningPoints(turningPoints==0) = []; %Deletes zero value
    
    %Do actual proper time evolution
    for i=2:POSITION_ARRAY_SIZE
        pos_t(i) = pos_t(i-1) + d_tau*(P_E/(1-2*BH_M/pos_r(i-1)));
        pos_phi(i) = pos_phi(i-1) + d_tau*(P_L/(pos_r(i-1)^2));
        %Check if particle is ingoing or outgoing and advance properly
        if ingoing_flag == 1
            pos_r(i) = pos_r(i-1) - d_tau*...
                (sqrt(P_E^2-(1-2*BH_M/pos_r(i-1))*(1+P_L^2/pos_r(i-1)^2)));
        else
            pos_r(i) = pos_r(i-1) + d_tau*...
                (sqrt(P_E^2-(1-2*BH_M/pos_r(i-1))*(1+P_L^2/pos_r(i-1)^2)));
        end
        
        %Check if particle crosses turning point & if so maintain current
        %radius and switch ingoing/outgoing flag
        for elm = turningPoints
            if (pos_r(i) - elm)*(pos_r(i-1) - elm) < 0
                pos_r(i) = pos_r(i-1);
                if ingoing_flag == 1
                    ingoing_flag = 0;
                    break
                else
                    ingoing_flag = 1;
                    break
                end
            end
        end
        
        %If particle is in Event horizon, mess stuff up
        if pos_r(i) < 2*BH_M
            pos_r(i:POSITION_ARRAY_SIZE) = ...
                2*BH_M*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_phi(i:POSITION_ARRAY_SIZE) = ...
                2*pi*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_t(i:POSITION_ARRAY_SIZE) = 0;
            return
        end
    end %End proper time progression loop

    %End of Schwarzchild treatment
elseif BH_Type == 0
    
    %Calculate initial coordinate velocities
    phi_dot = P_L/P_R^2;
    r_dot = sqrt(2*P_E-P_L^2/P_R^2+2*BH_M/P_R);
    if imag(r_dot) ~= 0
        fprintf('Invalid Energy/Angular Momentum Selection');
        pos_t = [];
        return;
    end
    
    r_double_dot = -P_L^2/P_R^3 + BH_M/P_R^2;
    
    %Circular orbit, muy bueno
    if r_double_dot == 0
        for i=2:POSITION_ARRAY_SIZE
            pos_t(i) = pos_t(i-1)+d_tau;
            pos_phi(i) = pos_phi(i-1) + d_tau*(P_L/P_R^2);
            pos_r(i) = P_R;
        end
        return
    end %end circular orbit
    
    
    %Solve for eccentricity and true anomaly in accordance with Wikipedia
    %Article on Kepler Orbits, wrt Determintation from Initial conditions
    p = (P_R*phi_dot)^2/BH_M;
    F = @(y) [y(1)*cos(y(2))+1-phi_dot/sqrt(BH_M/p) ...
        y(1)*sin(y(2))-r_dot/sqrt(BH_M/p)];
    
    eAndTheta = fsolve(F,[1 pi]);
    
    
    if eAndTheta(1) < 1 %Elliptic Orbit
        for i=2:POSITION_ARRAY_SIZE
            pos_t(i) = pos_t(i-1)+d_tau;
            pos_phi(i) = pos_phi(i-1) + d_tau*(P_L/pos_r(i-1)^2);
            pos_r(i) = p./(1+eAndTheta(1).*cos(pos_phi(i)-eAndTheta(2)));
        end
        return
        
    else %Parabolic or hyperbolic orbit
        
        turningPointS = roots([2*P_E,2*BH_M,-P_L^2]);
        turningPointS(imag(turningPointS) ~= 0) = [];
        turningPoint = max(turningPointS);

        for i=2:POSITION_ARRAY_SIZE
            pos_t(i) = pos_t(i-1)+d_tau;

            if ingoing_flag == 1
                pos_r(i) = pos_r(i-1) - d_tau *... 
                    sqrt(2*P_E-P_L^2/pos_r(i-1)^2+2*BH_M/pos_r(i-1));
                %pos_phi(i) = eAndTheta(2) + acos( 1/eAndTheta(1)* (p./pos_r(i) - 1));
                %If less than turning point, reverse
                if pos_r(i) < max(0,turningPoint)
                    pos_r(i) = pos_r(i-1);
                    ingoing_flag = 0;
                end

            else %ingoing_flag == 0
                pos_r(i) = pos_r(i-1) + d_tau *... 
                    sqrt(2*P_E-P_L^2/pos_r(i-1)^2+2*BH_M/pos_r(i-1));
                %pos_phi(i) = eAndTheta(2) - acos( 1/eAndTheta(1)* (p./pos_r(i) - 1));
            end

        end %Done handling pos_r/t values
        
        for i=2:POSITION_ARRAY_SIZE
            pos_phi(i) = pos_phi(i-1) + abs(...
                acos(1/eAndTheta(1)* (p./pos_r(i) - 1)) - ...
                acos( 1/eAndTheta(1)* (p./pos_r(i-1) - 1) ));
        end
        
        %pos_phi = eAndTheta(2) + acos( 1/eAndTheta(1)* (p./pos_r - 1) );
        %pos_phi = acos( 1/eAndTheta(1)* (p./pos_r - 1) );
            
    end % end of elliptic orbit handling
    return
    
elseif BH_Type == 2 %Kerr
    
    %Check for valid energy/Angular momentum values
    if (-1*K*BH_M/P_R + P_L^2/(2*P_R^2) + 0.5*(K-P_E^2)*(1+ (BH_L/P_R)^2) -...
            BH_M/P_R^3*(P_L-BH_L*P_E)^2) > 0
        fprintf('Invalid Energy/Angular Momentum Selection');
        pos_t = [];
        return;
    end
    
    %Create functions for E/L inversion wrt t_dot/phi_dot
    A = @(y) 1-2*BH_M./y;
    B = @(y) (y + BH_L^2./y).^2 - (y.^2+BH_L^2-2*BH_M.*y)*BH_L^2./y.^(-2);
    C = @(y) 2*BH_M*BH_L./y;
    
    t_dot = @(y) 1./(A(y).*B(y) + C(y).^2).* ...
        (P_E*B(y)-P_L*C(y));
    phi_dot = @(y) 1./(A(y).*B(y) + C(y).^2).* ...
        (P_E*C(y)+P_L*A(y));
    
    %Find Turning points based on tailored potential
    potentialPolynom = [0.5*(K-P_E^2), -1*K*BH_M, 0.5*(P_L^2+(K-P_E)^2*BH_L^2), ...
        -1*BH_M*(P_L-BH_L*P_E)^2];
    turningRoots = roots(potentialPolynom);
    turningPoints = zeros(1,3);
    %potentialPolynom = [(P_E^2-1) 2*BH_M -P_L^2 2*BH_M*P_L^2];
    
    
    for i = 1:length(turningRoots)
        if turningRoots(i) > 0 && imag(turningRoots(i)) == 0
            turningPoints(i) = turningRoots(i);
        end
    end
    turningPoints(turningPoints==0) = []; %Deletes zero value
    
    
    
    %Do actual proper time evolution
    for i=2:POSITION_ARRAY_SIZE
        pos_t(i) = pos_t(i-1) + d_tau*t_dot(pos_r(i-1));
        pos_phi(i) = pos_phi(i-1) + d_tau*phi_dot(pos_r(i-1));
        %Check if particle is ingoing or outgoing and advance properly
        if ingoing_flag == 1
            pos_r(i) = pos_r(i-1) - d_tau*...
                sqrt(-2*(-1*K*BH_M/pos_r(i-1) + P_L^2/(2*pos_r(i-1)^2) + ...
                0.5*(K-P_E^2)*(1 + (BH_L^2/pos_r(i-1)^2)) - BH_M/pos_r(i-1)^3*(P_L-BH_L*P_E)^2));
            if imag(pos_r(i)) ~= 0
                pos_r(i) = pos_r(i-2);
                ingoing_flag = 0;
            end
            
        else %outgoing
            pos_r(i) = pos_r(i-1) + d_tau*...
                sqrt(-2*(-1*K*BH_M/pos_r(i-1) + P_L^2/(2*pos_r(i-1)^2) + ...
                0.5*(K-P_E^2)*(1 + (BH_L^2/pos_r(i-1)^2)) - BH_M/pos_r(i-1)^3*(P_L-BH_L*P_E)^2));
            if imag(pos_r(i)) ~= 0
                pos_r(i) = pos_r(i-2);
                ingoing_flag = 1;
            end
        end
        
        %Check if particle crosses turning point & if so maintain current
        %radius and switch ingoing/outgoing flag
        for elm = turningPoints
            if (pos_r(i) - elm)*(pos_r(i-1) - elm) < 0
                pos_r(i) = pos_r(i-1);
                if ingoing_flag == 1
                    ingoing_flag = 0;
                    break
                else
                    ingoing_flag = 1;
                    break
                end
            end
        end
        
        %If particle is in Event horizon, mess stuff up
        if pos_r(i) < BH_M+sqrt(BH_M^2-(BH_L/BH_M)^2)
            pos_r(i:POSITION_ARRAY_SIZE) = ...
                2*BH_M*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_phi(i:POSITION_ARRAY_SIZE) = ...
                2*pi*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_t(i:POSITION_ARRAY_SIZE) = 0;
            return
        end
        
        %If particle is in Event horizon, mess stuff up
        %{
        if imag(pos_r(i)) ~= 0 && ingoing_flag ==0
            pos_r(i:POSITION_ARRAY_SIZE) = ...
                2*BH_M*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_phi(i:POSITION_ARRAY_SIZE) = ...
                2*pi*rand(1,POSITION_ARRAY_SIZE-i+1);
            pos_t(i:POSITION_ARRAY_SIZE) = 0;
            return
        end
        %}
        
    end %End proper time progression loop

    %End of Kerr treatment
    
end
