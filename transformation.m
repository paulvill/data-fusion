function [ T, Tu, Tl ] = transformation( indl, indu )
%TRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here

n = length(indl)+length(indu);
indtot = [indl,indu];
T = zeros(length(indu)+length(indl));
Tl = zeros(length(indl),length(indu)+length(indl));
Tu = zeros(length(indu),length(indu)+length(indl));

for i = 1:n,
    for j = 1:n,
        if j == indtot(i),
            T(i,j) = 1;
            if i <= length(indl),
                Tl(i,j) = 1;
            else
                Tu(i - length(indl),j) = 1;
            end
        end
    end
end

end

