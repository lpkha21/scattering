function f = fij(di,dj,ti,tj)
    yi = di * sin(ti);
    xi = di * cos(ti);
    yj = dj * sin(tj);
    xj = dj * cos(tj);
    
    f = mod(atan((yj - yi)/(xj - xi))+2*pi,2*pi);
end