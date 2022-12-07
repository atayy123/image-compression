function [yll, yhl, ylh, yhh] = fwt2d(x, h)
    yl1 = zeros(length(x), length(x)/2);
    yh1 = zeros(length(x), length(x)/2);

    yll = zeros(length(x)/2, length(x)/2);
    yhl = zeros(length(x)/2, length(x)/2);
    ylh = zeros(length(x)/2, length(x)/2);
    yhh = zeros(length(x)/2, length(x)/2);

    for i = 1:length(x)
        y1 = fwt(x(i,:), h);
        yl1(i,:) = y1(1:length(y1)/2);
        yh1(i,:) = y1(length(y1)/2+1:end);
    end

    for j = 1:size(yh1,2)
        yl2 = fwt(yl1(:,j), h);
        yh2 = fwt(yh1(:,j), h);
        
        yll(:,j) = yl2(1:length(yl2)/2);
        ylh(:,j) = yl2(length(yl2)/2+1:end);
        yhl(:,j) = yh2(1:length(yh2)/2);
        yhh(:,j) = yh2(length(yh2)/2+1:end);
    end

end