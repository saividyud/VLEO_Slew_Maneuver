function tauBody = rigid_body_total_torque(omegaBody, alphaBody, inertiaBody)
    tauBody = inertiaBody * alphaBody + cross(omegaBody, inertiaBody * omegaBody);
end
