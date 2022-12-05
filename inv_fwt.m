% inverse FWT
function [x] = inv_fwt(y, h)
    % find wavelet vector from scaling vector
    N = length(h);
    hw = zeros(N,1);
    for i = 0:N-1
        hw(i+1) = (-1)^i * h(N-i);
    end
    %disp(hw)
    % divide the high and low band
    y0 = y(1:length(y)/2);
    y1 = y(length(y)/2+1:end);

    % upsample the signals
    y0 = upsample(y0, 2);
    y1 = upsample(y1, 2);

    % extend the input signals
    y0_ext = [y0; y0; y0];
    y1_ext = [y1; y1; y1];

    % apply the filters
    x0 = conv(y0_ext, h, 'same');
    x0 = x0(length(y)+1:end-length(y));

    x1 = conv(y1_ext, hw, 'same');
    x1 = x1(length(y)+1:end-length(y));

    % merge the outputs
    x = x0 + x1;
    x = circshift(x,1);
end