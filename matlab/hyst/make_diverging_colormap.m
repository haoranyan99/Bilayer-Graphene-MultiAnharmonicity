function cmap = make_diverging_colormap()
% blue -> white -> red

    N = 256;
    half = floor(N/2);

    % blue to white
    blue = [0 0 0.8];
    white = [1 1 1];
    red = [0.8 0 0];

    cmap1 = [linspace(blue(1),white(1),half)', ...
             linspace(blue(2),white(2),half)', ...
             linspace(blue(3),white(3),half)'];

    cmap2 = [linspace(white(1),red(1),half)', ...
             linspace(white(2),red(2),half)', ...
             linspace(white(3),red(3),half)'];

    cmap = [cmap1; cmap2];
end