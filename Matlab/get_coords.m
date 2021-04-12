function coords = get_coords(coord_dirs, res, im0)

coords = zeros(4);
coords(4,4) = 1.0;
for i = 1:3
    j = coord_dirs(i);
    coords(i,j) = res(j);
    coords(i,4) = -size(im0,j) * res(j) / 2.0;
end

%Fix matlab 1 indexing
Q = zeros(4);
Q(1:3,4) = ones(3,1);
coords = inv(inv(coords)+Q);

end

