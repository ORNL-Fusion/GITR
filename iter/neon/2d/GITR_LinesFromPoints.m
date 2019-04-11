function lines = GITR_LinesFromPoints(r,z, openClosed)

switch openClosed
    case 'closed'
        nPoints = length(r);
        lines = zeros(nPoints,7);
        lines(:,1) = r;
        lines(:,2) = z;
        lines(1:end-1,3) = r(2:end);
        lines(1:end-1,4) = z(2:end);
        lines(end,3) = r(1);
        lines(end,4) = z(1);
    case 'open'
        nPoints = length(r)-1;
        lines = zeros(nPoints,7);
        lines(:,1) = r(1:(nPoints));
        lines(:,2) = z(1:(nPoints));
        lines(1:(nPoints-1),3) = lines(2:end,1);
        lines(1:(nPoints-1),4) = lines(2:end,2);
        lines(end,3) = r(end);
        lines(end,4) = z(end);
end

tol = 1e12;
tol_small = 1e-12;
for i=1:nPoints
    if (lines(i,4) - lines(i,2)) == 0
        lines(i,5) = 0;
        lines(i,6) = lines(i,2);
    else if (lines(i,3) - lines(i,1)) == 0;
            lines(i,5) = sign(lines(i,4) - lines(i,2))*tol;
            lines(i,6) = tol;
        else
        lines(i,5) = (lines(i,4) - lines(i,2))/(lines(i,3) - lines(i,1));
        lines(i,6) = -lines(i,5)*lines(i,1) + lines(i,2);
        end
    end
end

lines(:,7) = sqrt((lines(:,3) - lines(:,1)).^2 + (lines(:,4) - lines(:,2)).^2);