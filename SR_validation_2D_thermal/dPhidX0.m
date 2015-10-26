function out = dPhidX0(x, x0, h)
r = (x - x0) / h;
out = zeros(size(r,1), size(r,2));

%% 2-Point DELTA Function
for i=1:size(r, 1)
    for j=1:size(r, 2)
        if (r(i, j) < 1) && (r(i, j) > 0)
            out(i, j) = -1;
        elseif (r(i, j) < 0) && (r(i, j) > -1)
            out(i, j) = 1;
        end
    end
end
out = -out / h;

%% 4-Point DELTA Function
% for i=1:size(r, 1)
%     for j=1:size(r, 2)
%         if abs(r(i, j))<= 1
%             out(i, j) = (3 - 2 * abs(r(i, j)) + sqrt(1 + 4 * abs(r(i, j)) - 4 * r(i, j)^2)) / 8;
%         elseif (abs(r(i, j)) <= 2) && (abs(r(i, j)) >= 1)
%             out(i,j) = (5 - 2 * abs(r(i, j)) - sqrt(-7 + 12 * abs(r(i, j)) - 4 * r(i,j)^2)) / 8;
%         end
%     end
% end