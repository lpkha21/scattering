function f = fij(di,dj,ti,tj)
    yi = di * sin(ti);
    xi = di * cos(ti);
    yj = dj * sin(tj);
    xj = dj * cos(tj);
    
    f = vpa(atan2(yj-yi,xj-xi),500);
    
    if f < 0
       f = vpa(f + 2*pi,500);
    end
    
    
end