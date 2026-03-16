function Xd = sat_dynamics_openloop(t, X, params, tauControlInterp, tauAeroInterp)
    qNorm = norm(X(7:10));
    if qNorm < 1e-10 || ~isfinite(qNorm)
        Xd = zeros(13, 1);
        return;
    end
    X(7:10) = X(7:10) / qNorm;
    tauBody = reshape(tauControlInterp(t), [], 1) + reshape(tauAeroInterp(t), [], 1);
    Xd = vleo.dynamics.sat_dynamics_base(t, X, params.mu, [], tauBody, params.I_CB, false);
end
