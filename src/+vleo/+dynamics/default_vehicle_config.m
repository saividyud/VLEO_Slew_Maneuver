function vehicle = default_vehicle_config()
    persistent cachedVehicle

    if isempty(cachedVehicle)
        cachedVehicle = struct();
        cachedVehicle.objFilePath = vleo.util.geometry_asset_path('6U CubeSat.obj');
        cachedVehicle.stlFilePath = vleo.util.geometry_asset_path('6U CubeSat.STL');
        cachedVehicle.mass = 83;
        cachedVehicle.inertiaBody = 2 / 5 * cachedVehicle.mass * (0.58 / 2)^2 * eye(3);
    end

    vehicle = cachedVehicle;
end
