function [ cccoeff ] = CCC(x,y)
    m1 = mean(x);
    m2 = mean(y);
    s1 = var(x)*(length(x)-1)/length(x);
    s2 = var(y)*(length(y)-1)/length(y);
    s12 = sum((x-m1).*(y-m2))/length(x);
    cccoeff = 2*s12 / (s1+s2+(m1-m2).^2);
end

