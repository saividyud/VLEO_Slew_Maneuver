function vehicle = default_vehicle_config()
    persistent cachedVehicle

    if isempty(cachedVehicle)
        cachedVehicle = struct();
        cachedVehicle.objFilePath = vleo.util.geometry_asset_path('6U CubeSat.obj');
        cachedVehicle.stlFilePath = vleo.util.geometry_asset_path('6U CubeSat.STL');
        cachedVehicle.mass = 12;
        c = .3;
        b = .2; 
        a = .1;
        cachedVehicle.inertiaBody = 1 / 12 * cachedVehicle.mass * [ (b^2 + c^2) 0 0 ; 0 (a^2 + c^2) 0 ; 0 0 (a^2 + b ^2)];
    end

    vehicle = cachedVehicle;
end
