function Xd = Sat_template(t, X, varargin)
    % Extract control torque from varargin
    % Handles both: Sat_template(t, X, tau) and Sat_template(t, X, tau, params, 'useJ2', false, ...)
    
    if nargin >= 3
        if isstruct(varargin{1})
            % First arg is control torque, second is params
            tau_body = varargin{1};
            params = varargin{2};
        else
            tau_body = varargin{1};
            % Extract params from varargin or use defaults
            params.I_CB = 2/5 * 83.6 * 0.29^2 * eye(3);
            params.mu = 3.986004418e14;
        end
    else
        error('Insufficient arguments');
    end
    
    % Rest of dynamics computation...
    Xd = zeros(13, 1);
    r_eci = X(1:3);
    v_eci = X(4:6);
    q = X(7:10);
    omega_eci = X(11:13);
    q = q / norm(q);
    
    Xd(1:3) = v_eci;
    r_mag = norm(r_eci);
    mu = params.mu;
    a_grav = -mu * r_eci / r_mag^3;
    Xd(4:6) = a_grav;
    
    R_b2i = quat_to_dcm(q);
    R_i2b = R_b2i';
    omega_body = R_i2b * omega_eci;
    
    Omega_matrix = [ 0, omega_body(3), -omega_body(2), omega_body(1);
                   -omega_body(3), 0, omega_body(1), omega_body(2);
                    omega_body(2), -omega_body(1), 0, omega_body(3);
                   -omega_body(1), -omega_body(2), -omega_body(3), 0];
    Xd(7:10) = 0.5 * Omega_matrix * q;
    
    I_body = params.I_CB;
    Iomega = I_body * omega_body;
    omega_cross_Iomega = cross(omega_body, Iomega);
    alpha_body = I_body \ (tau_body - omega_cross_Iomega);
    alpha_eci = R_b2i * alpha_body;
    Xd(11:13) = alpha_eci;
end

function R = quat_to_dcm(q)
    qx = q(1); qy = q(2); qz = q(3); qw = q(4);
    R = [1 - 2*(qy^2 + qz^2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy);
         2*(qx*qy + qw*qz), 1 - 2*(qx^2 + qz^2), 2*(qy*qz - qw*qx);
         2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx^2 + qy^2)];
end
