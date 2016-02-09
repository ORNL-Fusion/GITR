function b2_cplot(crx,cry,data,datamin);

if nargin < 4 
    datamin = -1e100;
end

s=size(crx);
figure
for i = 1:s(1)
    for j = 1:s(2)
        if data(i,j) > datamin
            patch([crx(i,j,1) crx(i,j,2) crx(i,j,4) crx(i,j,3)],[cry(i,j,1) cry(i,j,2) cry(i,j,4) cry(i,j,3)], ...
                data(i,j),'edgecolor','none')
        end
    end
end
axis equal
