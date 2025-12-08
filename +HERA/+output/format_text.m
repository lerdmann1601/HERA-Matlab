function output = format_text(text, width, alignment)
% FORMAT_TEXT - Helper function to format text with padding and alignment.
%
% Syntax:
%   output = format_text(text, width, alignment)
%
% Description:
%   Formats a string to a specific width with the specified alignment ('l', 'c', 'r').
%   Used for generating aligned table outputs in the console.
%
% Inputs:
%   text      - The string to format.
%   width     - The total width of the output string (including padding).
%   alignment - 'l' (left), 'c' (center), or 'r' (right).
%
% Outputs:
%   output    - The formatted string with requested padding.
%
% Author: Lukas von Erdmannsdorff

    % Width without the buffer padding of 2, as this is only relevant for 'c'
    text_len = strlength(text);
    
    switch alignment
        case 'r'
            % Right-aligned: Uses the total width 'width' minus the text length.
            % The 2 padding characters are included here but displayed as spaces to the left of the text.
            output = [repmat(' ', 1, width - text_len), text];
        case 'l'
            % Left-aligned: Uses the total width 'width' minus the text length.
            % The 2 padding characters are included here, displayed as spaces to the right of the text.
            output = [text, repmat(' ', 1, width - text_len)];
        otherwise % 'c' (Centered)
            % Centered: Recalculates padding to center the text.
            padding = width - text_len;
            padding_left = floor(padding / 2);
            padding_right = ceil(padding / 2);
            output = [repmat(' ', 1, padding_left), text, repmat(' ', 1, padding_right)];
    end
end
