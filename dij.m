function d = dij(di,dj,ti,tj)
    d = sqrt(di*di + dj*dj - 2*di*dj*cos(ti-tj));
end