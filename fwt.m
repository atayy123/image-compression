% FWT filter bank with direct implementation
function [y] = fwt(x, h)

    % find wavelet vector from scaling vector
    N = length(h);
    hw = zeros(N,1);
    for i = 1:N
        hw(i) = (-1)^i * h(N-i+1);
    end

    %disp(hw)
    % FWT direct implementation
    
    % extend the input signal
    x_ext = [x; x; x];

    % apply the filters
    y0 = conv(x_ext, h, "same");
    y0 = y0(length(x)+1:end-length(x));
    y0 = downsample(y0, 2);

    y1 = conv(x_ext, hw, 'same');
    y1 = y1(length(x)+1:end-length(x));
    y1 = downsample(y1, 2);

    y = [y0; y1];
end