% Approximate comparison function
function result = isapprox(a, b, atol)
    result = all(abs(a - b) < atol);
end