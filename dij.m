function d = dij(di,dj,ti,tj)
    yi = di * sin(ti);
    xi = di * cos(ti);
    yj = dj * sin(tj);
    xj = dj * cos(tj);
    
    d = sqrt((xi-xj)^2 + (yi-yj)^2);
end