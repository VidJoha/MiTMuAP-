function Z = blondel_slicnost(A, B, Z0, maxIter, tol)
  Z = Z0 / norm(Z0, 'fro');
    for it = 1:maxIter
        Zold = Z;
        for step = 1:2
            Z = B*Z*A' + B'*Z*A;
            Z = Z / norm(Z, 'fro');
        end
        if norm(Z - Zold, 'fro') < tol
            break;
        end
    end
end
