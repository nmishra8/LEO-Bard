function [ error ] = accuracy_rss( x,y )
    residual = x-y;
    n = length(x);
    rss = sum(residual.*residual); 
    tss = sum((x-mean(x)).*(x-mean(x)));
    residualsquare = 1-rss/tss;
    error = 1-(1-residualsquare)*(n-1)/(n-2);
end

