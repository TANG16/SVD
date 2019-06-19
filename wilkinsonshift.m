function mu = wilkinsonshift( a1, b, a2 )

d = ( a1 - a2 ) / 2; % flops: 2

if d == 0

    if a2 > 0
        mu = a2 + abs( b );
        else
        mu = a2 - abs( b );
    end
else

    mu = a2 - b^2 / ( d + sign(d)*sqrt( d^2 + b^2 ) );
end

end