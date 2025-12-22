function [X_analytical] = Analytical_Propagator(X0, t, params)
% Analytical_Propagator propagates satellite orbit using Kepler's equations
%
% This function uses analytical solutions to the two-body problem to
% propagate orbital state. It does NOT include perturbations (J2, drag, etc.)
% and is used as a benchmark for comparison with numerical propagators.
%
% Parameters
% ----------
% X0 : 6x1 or 13x1 vector
%     Initial state vector containing:
%     - Position in ECI (1:3) [m]
%     - Velocity in ECI (4:6) [m/s]
%     If 13x1, only first 6 elements are used
% t : float or array
%     Time or array of times to propagate to [s]
% params : struct
%     Physical parameters with field:
%     - mu: Gravitational parameter [m^3/s^2]
%
% Returns
% -------
% X_analytical : 6xN array
%     State vectors at requested times (position + velocity only)
%     Each column is [x; y; z; vx; vy; vz] at corresponding time
%
% Algorithm
% ---------
% Uses Kepler's universal variable formulation with Newton-Raphson iteration
% to solve the time-of-flight equation. This is valid for all orbit types
% (elliptic, parabolic, hyperbolic).
%
% References
% ----------
% Vallado, D. A. "Fundamentals of Astrodynamics and Applications"

    % Input validation
    assert(length(X0) >= 6, 'Initial state must have at least 6 elements');
    assert(isfield(params, 'mu'), 'params must contain mu');

    % Extract orbital elements from initial state
    r0_vec = X0(1:3);
    v0_vec = X0(4:6);
    mu = params.mu;

    % Initialize output
    n_times = length(t);
    X_analytical = zeros(6, n_times);

    % Calculate initial orbit parameters
    r0 = norm(r0_vec);
    v0 = norm(v0_vec);

    % Radial velocity
    vr0 = dot(r0_vec, v0_vec) / r0;

    % Specific orbital energy
    epsilon = v0^2/2 - mu/r0;

    % Semi-major axis
    if abs(epsilon) > 1e-10
        a = -mu / (2*epsilon);
    else
        a = inf;  % Parabolic orbit
    end

    % Propagate for each time
    for i = 1:n_times
        dt = t(i);

        if abs(dt) < 1e-10
            % No time elapsed
            X_analytical(:, i) = [r0_vec; v0_vec];
        else
            % Solve universal Kepler's equation using Newton-Raphson
            [r_vec, v_vec] = kepler_universal(r0_vec, v0_vec, dt, mu);
            X_analytical(:, i) = [r_vec; v_vec];
        end
    end

end


function [r_vec, v_vec] = kepler_universal(r0_vec, v0_vec, dt, mu)
% kepler_universal solves Kepler's problem using universal variables
%
% This function uses the f and g functions with universal anomaly

    r0 = norm(r0_vec);
    v0 = norm(v0_vec);

    % Radial velocity
    vr0 = dot(r0_vec, v0_vec) / r0;

    % Reciprocal of semi-major axis
    alpha = 2/r0 - v0^2/mu;

    % Initial guess for universal anomaly
    if alpha > 1e-6
        % Elliptical orbit
        X = sqrt(mu) * dt * alpha;
    else
        % Parabolic or hyperbolic
        X = 0;
    end

    % Newton-Raphson iteration
    max_iter = 100;
    tol = 1e-10;

    for iter = 1:max_iter
        X2 = X^2;
        psi = alpha * X2;

        % Stumpff functions
        [C, S] = stumpff(psi);

        % Time equation
        F = r0*vr0/sqrt(mu) * X2 * C + (1 - alpha*r0)*X^3*S + r0*X - sqrt(mu)*dt;

        % Derivative
        dFdX = r0*vr0/sqrt(mu) * X * (1 - alpha*X2*S) + ...
               (1 - alpha*r0)*X2*C + r0;

        % Newton-Raphson update
        X_new = X - F/dFdX;

        % Check convergence
        if abs(X_new - X) < tol
            X = X_new;
            break;
        end

        X = X_new;
    end

    % Calculate f and g functions
    X2 = X^2;
    psi = alpha * X2;
    [C, S] = stumpff(psi);

    f = 1 - X2/r0 * C;
    g = dt - X^3/sqrt(mu) * S;

    % Position at time t
    r_vec = f*r0_vec + g*v0_vec;
    r = norm(r_vec);

    % f_dot and g_dot
    f_dot = sqrt(mu)/(r*r0) * X * (alpha*X2*S - 1);
    g_dot = 1 - X2/r * C;

    % Velocity at time t
    v_vec = f_dot*r0_vec + g_dot*v0_vec;

end


function [C, S] = stumpff(psi)
% stumpff calculates Stumpff functions C(psi) and S(psi)
%
% These functions are used in the universal variable formulation

    if psi > 1e-6
        % Elliptical orbit
        sqrt_psi = sqrt(psi);
        C = (1 - cos(sqrt_psi)) / psi;
        S = (sqrt_psi - sin(sqrt_psi)) / sqrt_psi^3;
    elseif psi < -1e-6
        % Hyperbolic orbit
        sqrt_neg_psi = sqrt(-psi);
        C = (1 - cosh(sqrt_neg_psi)) / psi;
        S = (sinh(sqrt_neg_psi) - sqrt_neg_psi) / (-psi)^(3/2);
    else
        % Parabolic or near-parabolic orbit (use series expansion)
        C = 1/2 - psi/24 + psi^2/720;
        S = 1/6 - psi/120 + psi^2/5040;
    end

end
