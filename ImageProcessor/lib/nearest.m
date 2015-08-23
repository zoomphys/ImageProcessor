function index = nearest(v,x)
% find the index in the vector v (or row index if v is a matrix) that contains value (or row of values) that is closes to x

% make v a column vector if v is a row vector and x is a scalar 
if isrow(v) && isscalar(x)
    v=v';
end

% replicate x to match the size of v and take the difference with v
diff = v-repmat(x,length(v(:,1)),1);

% take the sum along each row
distance = sum(diff.^2,2);

[value,index] = min(distance);

