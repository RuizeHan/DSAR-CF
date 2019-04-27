function Y = normalizing(X,a,b)
    minX=min(min(X));
    maxX=max(max(X));
    Y=(X-minX)*((b-a)/(maxX-minX))+a;
end