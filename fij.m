function f = fij(di,dj,ti,tj)
    yi = di * sin(ti);
    xi = di * cos(ti);
    yj = dj * sin(tj);
    xj = dj * cos(tj);
    
    f = atan2(yj-yi,xj-xi);
    
    if f < 0
       f = f + 2*pi;
    end
end