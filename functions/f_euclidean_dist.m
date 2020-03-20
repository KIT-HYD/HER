function [d] = f_euclidean_dist(x1, y1, x2, y2)
% Calculate the euclidean distance between the coords

%x1,x2 = coordenates X of the points (it can be a vector)
%y1,y2 = coordenates Y of the points (it can be a vector)
d = (sqrt((x1-x2).^2 + (y1-y2).^2));

end

