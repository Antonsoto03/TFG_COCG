function [x,iter,resvec] = cocg(A,b,x0,tol,maxIter)
%------------------------------------------------------------------

    x      = x0;
    r      = b - A*x;                 % r_0
    p      = r;                       % p_0
    rho    = conj(r)'*r;             % ⟨r0, r̄0⟩
    res0   = norm(b);  if res0==0, res0 = 1; end
    resvec = zeros(maxIter+1,1);      % guarda residuales
    resvec(1) = norm(r);

    for k = 0:maxIter-1
        q    = A*p;            % A p̄_k
        denom = conj(q)'*p;               % ⟨p_k , A p̄_k⟩   (p.' : transpuesta SÍMPL.)
        if abs(denom)<tol
            warning('COCG:breakdown','p^T A p ≈ 0 (k = %d).',k); break
        end
        alpha = rho / denom;          % α_k

        %----- x_{k+1} , r_{k+1}  -----------------------------------
        x = x + alpha * p;            % x_{k+1}
        r = r - alpha * q;           % r_{k+1}

        res = norm(r);
        resvec(k+2) = res;
        if res <= tol*res0            % convergencia
            iter   = k+1;
            resvec = resvec(1:iter+1);
            return
        end

       rho_new = conj(r)'*r;        % ⟨r_{k+1}, r̄_{k+1}⟩
        if k==0
            % para k = 0, el denominador usa r_{-1}; definimos β_0 = 0
            beta = 0;
        else
            beta = rho_new / rho;     % β_k
        end

        %----- p_{k+1}  ---------------------------------------------
        p   = r + beta * p;           % p_{k+1}
        rho = rho_new;                % listo para la siguiente vuelta
    end

    iter   = k;
    resvec = resvec(1:iter+1);
    warning('COCG:NoConverge','No converge tras %d iteraciones',iter);
end
