%Code for homework #4 in Orbital Mechanics
%Solve Lamberts problem when only given dt,r1, and r2
%Outputs should be v1 and v2 to correspond to r1 and r2


%% Test Cases 
%% Case 1
% Variables
r0 = [0.5 0.6 0.7];
r1 = [0.0 1.0 0.0];
dt = 0.9667663;
u = 1.0;
dir = -1;
% Answer Print
[v0,v1] = Lambert(r0, r1, dt, u, dir);
fprintf('v0 = [');
fprintf('%g ', v0);
fprintf(']\n');
fprintf('v0 actual = [ -0.631 -1.114 -0.883]\n'); 
fprintf('v1 = [');
fprintf('%g ', v1);
fprintf(']\n');
fprintf('v1 actual = [0.1787 1.554 0.250]\n\n');
%% Case 2
% Variables
r02 = [1 0 0];
r12 = [1 0.125 0.125];
dt2 = 0.125;
u = 1.0;
dir2 = 1;
% Answer Print
[v0,v1] = Lambert(r02, r12, dt2, u, dir2);
fprintf('v0 = [');
fprintf('%g ', v0);
fprintf(']\n');
fprintf('v0 actual = [0.0619 1.003 1.003]\n'); 
fprintf('v1 = [');
fprintf('%g ', v1);
fprintf(']\n');
fprintf('v1 actual = [-0.0609 0.995 0.995]\n');








%% Answer functions
function [v0, v1] = Lambert(r0, r1, dt, u, dir)
    %% Variable creation
R0 = r0;
R1 = r1;
R0mag = norm(R0);
R1mag = norm(R1);
Dt = dt;
mu = u;

    %% Simple Equations
dtheta = acos((dot(R0,R1))/R0mag/R1mag);
if dir == 1
    sign = 1;
else
    sign = -1;
    dtheta = 2*pi - dtheta;
end
A = sign*sqrt(R0mag*R1mag*(1+cos(dtheta)));

    %% Newton method to determine X
    z = NewtonMethod(Dt, A, mu, R0mag, R1mag);
    % take answer from Newton method and solve
        s_z = (sqrt(z)-sin(sqrt(z)))/(z^(3/2)); 
        c_z = (1 - cos(sqrt(z)))/z;
        y_z = R0mag + R1mag - A*(1-z*s_z)/(sqrt(c_z));
        X = sqrt(y_z/c_z);
        f = 1 - (X^2./R0mag).*c_z;
        g = Dt - (X^3*s_z)/sqrt(mu);
        v0 = (R1 - f.*R0) / g;
        df = (sqrt(mu)./(R0mag*R1mag)) *(X*(z*s_z-1));
        dg = (df * g + 1) / f;
        v1 = df.*R0 + dg.*v0;
end



%% Newton Method function
function [z] = NewtonMethod(dt, A, mu, r0, r1)
 max_iter = 100;
    tol = 1e-12;
    z = 0;
    for iter = 1:max_iter
        % calculate z using u(z) and v(z) equations
        if z == 0
            c_z = 1/2;
            s_z = 1/6;
            dc_z = -1/24;
            ds_z = -1/120;
        else
            s_z = (sqrt(z)-sin(sqrt(z)))/z^(3/2); 
            c_z = (1-cos(sqrt(z)))/z;
            ds_z = 1/(2*z)*(c_z-3*s_z);
            dc_z = 1/(2*z)*(1-z*s_z-2*c_z);
        end
            y_z = r0 + r1 - A*(1-z*s_z)/(sqrt(c_z));
            x = sqrt(y_z/c_z);
            u_z = 1/sqrt(mu)* ((x^3)*s_z + A*sqrt(y_z)) -dt;
            du_z = 1/sqrt(mu)* ((x^3)*(ds_z-(3*s_z*dc_z)/(2*c_z)) +(A/8* ((3*s_z*sqrt(y_z))/c_z +(A/x))));
        
        % Check for convergence
        if abs(u_z) < tol || iter == max_iter
            break;
        end
        % Update z to new value
        z =z-((u_z)/du_z);
    end
end