function tau = attitude_pd_controller(~, X)
    c = .3;
    b = .2; 
    a = .1;
    inertiaBody = 1 / 12 * 12 * [ (b^2 + c^2) 0 0 ; 0 (a^2 + c^2) 0 ; 0 0 (a^2 + b ^2)];
    qReference = [1; 0; 0; 0];
    omegaReference = [0; 0; 0];
    omegaDotReference = [0; 0; 0];
    dampingGain = 0.06 * eye(3);
    proportionalGain = 0.0014 * eye(3);

    deltaOmega = X(11:13) - omegaReference;
    tau = -proportionalGain * (X(8:10) - qReference(2:4));
    tau = tau - dampingGain * deltaOmega;
    tau = tau + inertiaBody * omegaDotReference;
    tau = tau - cross(X(11:13), omegaReference);
end
