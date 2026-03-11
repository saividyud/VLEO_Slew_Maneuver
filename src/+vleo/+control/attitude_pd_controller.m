function tau = attitude_pd_controller(~, X)
    inertiaBody = 2 / 5 * 83 * (0.58 / 2)^2 * eye(3);
    qReference = [1; 0; 0; 0];
    omegaReference = [0; 0; 0];
    omegaDotReference = [0; 0; 0];
    dampingGain = 0.1 * eye(3);
    proportionalGain = 0.002 * eye(3);

    deltaOmega = X(11:13) - omegaReference;
    tau = -proportionalGain * (X(8:10) - qReference(2:4));
    tau = tau - dampingGain * deltaOmega;
    tau = tau + inertiaBody * omegaDotReference;
    tau = tau - cross(X(11:13), omegaReference);
end
