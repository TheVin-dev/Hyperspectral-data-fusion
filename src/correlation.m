function CC = correlation(A,B)

Ahat = average(A);
Anorm = (A(:)-Ahat);
Bhat = average(B);
Bnorm = (B(:)-Bhat);
num = sum(Anorm .* Bnorm,'all');

den = sum(Anorm.^2) * sum(Bnorm.^2);

CC = num / sqrt(den);


end

function ave = average(x)
    ave = sum(x(:))/numel(x); 
end