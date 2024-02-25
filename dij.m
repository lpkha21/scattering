function d = dij(di,dj,ti,tj)
    xi = di * cos(ti);
    yi = di * sin(ti);
    xj = dj * cos(tj);
    yj = dj * sin(tj);
    
    d = vpa(sqrt((xi-xj)^2 + (yi-yj)^2),500);
end