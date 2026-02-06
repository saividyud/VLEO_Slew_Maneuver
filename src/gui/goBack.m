%% Helper Function: The Back Logic
function goBack(currentFig, mainFig)
    delete(currentFig);      % Destroy the sub-window
    mainFig.Visible = 'on';  % Reveal the main menu
end