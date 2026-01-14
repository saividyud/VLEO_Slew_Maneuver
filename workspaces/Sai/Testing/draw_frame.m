function arrows = draw_frame(vectors, color, parent)
    
    arrows = cell(1, 3);
    for i = 1:3
        arrows{i} = quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), color, 'LineWidth', 2, Parent=parent);
    end

end