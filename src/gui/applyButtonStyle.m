%% Helper Function: Styling
% We use this to apply Times New Roman to buttons to avoid repeating code
function applyButtonStyle(btn)
    lightGrayColor = [211, 211, 211] / 256;
    
    btn.FontName = 'Times New Roman';
    btn.FontSize = 14;
    btn.FontWeight = 'bold';
    btn.BackgroundColor = lightGrayColor;
end

