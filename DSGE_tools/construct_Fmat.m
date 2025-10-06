function F = construct_Fmat(nn,nl)

%nn = 5; % number of vars in Sig0
%nl = 4; % one plus number of lags (e.g., for Sig^2k, nl=2k+1)

% we want to construct F with dimensions (nl*nn)^2 x nl*nn^2, such that
% vec(Sig^nl) = F* [vec(Sig0); vec(Sig1); ...]

%F = zeros(nl^2*nn^2,nl*nn^2);
F  = spalloc(nl^2*nn^2,nl*nn^2,nn^2*nl);
for maincol = 1:nl
    
    for matcol = 1:nn
        
        for mainrow = 1:nl
            
            pos = (maincol-1)*nn^2*nl;
            pos = pos + (matcol-1)*nn*nl;
            pos = pos + (mainrow-1)*nn;
            
            if mainrow <= maincol
                
                F(pos+(1:nn),:) = ncolj(matcol, maincol-mainrow+1 , nn);
                
            else
                
                F(pos+(1:nn),:) = nrowj(matcol, mainrow-maincol+1 , nn);
                                
            end
            
        end
        
    end
    
    F = sparse(F);
end




%A = rand(nn); Avec=vec(A);

    function mat = ncolj(col,j,nn) % selects n-th column of Sig_j
        mat = zeros(nn,nl*nn^2);
        mat(:,nn^2*(j-1)+(1:nn^2))=ncol(col,nn);
    end

    function mat = nrowj(row,j,nn) % selects n-th row of Sig_j
        mat = zeros(nn,nl*nn^2);
        mat(:,nn^2*(j-1)+(1:nn^2))=nrow(row,nn);
    end


    function sel = ncol(col,nn) % selects n-th column of some Sig_j
        %sel = zeros(nn,nn^2);
        %sel(:,nn*(col-1)+(1:nn)) = eye(nn);
        sel = kron(evec(col,nn),eye(nn));
    end

     function sel = nrow(row,nn) % selects n-th row of some Sig_j
        sel = kron(eye(nn),evec(row,nn));
     end

    function e = evec(col,nn)
        e = zeros(1,nn);
        e(col) = 1;
    end
    
    
    
end