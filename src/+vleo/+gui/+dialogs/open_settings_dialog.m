function open_settings_dialog(subFigHandle, buttonType)
    switch buttonType
        case 'Set Simulation Parameters'
            vleo.gui.dialogs.set_simulation_parameters(subFigHandle);
        case 'Set Orbital Parameters'
            vleo.gui.dialogs.set_orbital_parameters(subFigHandle);
        case 'Set Attitude Parameters'
            vleo.gui.dialogs.set_attitude_parameters(subFigHandle);
        case 'Set Environmental Parameters'
            vleo.gui.dialogs.set_environmental_parameters(subFigHandle);
        case 'Set Controller Parameters'
            vleo.gui.dialogs.set_controller_parameters(subFigHandle);
        case 'Set Aero/Control Modes'
            vleo.gui.dialogs.set_aero_control_modes(subFigHandle);
        otherwise
            error('VLEO_Slew_Maneuver:UnknownButtonType', 'Unknown button type: %s', buttonType);
    end
end
