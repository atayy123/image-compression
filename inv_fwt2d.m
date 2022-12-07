function [x] = inv_fwt2d(yll, yhl, ylh, yhh, h)
    yl = zeros(2*length(yll),length(yll));
    yh = zeros(2*length(yll),length(yll));
 
    for i = 1:size(yll, 2)
        yl(:,i) = inv_fwt([yll(:,i); ylh(:,i)], h);
        yh(:,i) = inv_fwt([yhl(:,i); yhh(:,i)], h);
    end

    x = zeros(2*length(yll),2*length(yll));
    for j = 1:length(yl)
        x(j,:) = inv_fwt([yl(j,:) yh(j,:)], h);
    end
end