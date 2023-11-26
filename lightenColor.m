function lighterColor = lightenColor(color, amount)
    % color: RGB triplet [R, G, B]
    % amount: A value between 0 and 1, where 0 is the original color and 1 is white.

    if amount < 0 || amount > 1
        error('Amount must be between 0 and 1');
    end

    % Interpolate between the original color and white
    lighterColor = color + (1 - color) * amount;
end