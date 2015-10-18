function out = phi(x, x0, h)
r = (x - x0) / h;
out = zeros(size(r,1), size(r,2));
for i=1:size(r, 1)
    for j=1:size(r, 2)
        if abs(r(i, j))<= 1
            out(i, j) = 1 - abs(r(i, j));
        end
    end
end