function [tspan, xHistory] = propagate_bidirectional(odeFun, X0, tStart, tEnd, dt, odeOpts)
    tspanFwd = 0:dt:tEnd;
    [~, xFwd] = ode45(odeFun, tspanFwd, X0, odeOpts);

    tspanBwd = 0:-dt:tStart;
    [~, xBwd] = ode45(odeFun, tspanBwd, X0, odeOpts);

    xHistory = [flipud(xBwd(2:end, :)); xFwd];
    tspan = [fliplr(tspanBwd(2:end)), tspanFwd];
end
