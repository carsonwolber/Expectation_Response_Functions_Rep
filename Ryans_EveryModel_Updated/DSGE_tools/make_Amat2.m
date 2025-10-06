function [A1,A2,A3,A4,A] = make_Amat(K,bk_wghts,cross)

%Window of BK filter
k = (length(bk_wghts)-1)/2;


%These are the weights for level-detrended
A1 = zeros(K+1,length(bk_wghts)+K);
F = construct_Fmat(1,length(bk_wghts)+K);
for kk = 0:K
   A1(kk+1,:) = kron([zeros(1,kk),bk_wghts,zeros(1,K-kk)],[bk_wghts,zeros(1,K)])*F;
end


%Weights for the growing-detrended variables
J = 2*k+2*k+10;
S = ones(1,J+1);
Astar = zeros(J+1,length(bk_wghts)+J);
for jj = 1:2*k+1
    Astar(jj,jj:(jj+length(bk_wghts)-1)) = bk_wghts; 
end
Astar1 = Astar(1,:);  
F = construct_Fmat(1,length(bk_wghts)+J);
A2 = zeros(K+1,size(F,2));
A3 = zeros(K+1,size(F,2));
A4 = zeros(K+1,size(F,2));

%Formula in notes
SAstar_trunc = S*Astar;
SAstar_trunc(2*k+1:end) = 0;
SAstar_truncl = SAstar_trunc;
for kk = 1:K+1
    
    %Growth on growth weights
    A2(kk,:) = kron(SAstar_trunc,SAstar_truncl)*F;
    
    %Cross weights bottom left
    A3(kk,:) = kron(SAstar_trunc,Astar(1,:))*F;
    
    %Cross weights top right
    A4(kk,:) = kron(SAstar_truncl,Astar1 )*F;
    
    %Shifting timing
    Astar = [zeros(J+1,1),Astar(:,1:end-1)];
    SAstar_truncl = S*Astar;
    SAstar_truncl(kk+2*k+1:end) = 0;
end

cc = cross'*cross;
%Put in most covenience dimensions
A = zeros(K+1,2*k+K+1,length(cross),length(cross));
for jj = 1:K+1
    for kk = 1:(2*k+K+1)
        for ii =1:length(cross)
            for mm = 1:length(cross)
                
                if cc(ii,mm) ==4
                    A(jj,kk,ii,mm) = A2(jj,kk);
                elseif cc(ii,mm)==1
                    A(jj,kk,ii,mm) = A1(jj,kk);
                elseif cc(ii,mm)==2 && ii<mm
                    A(jj,kk,ii,mm) = A4(jj,kk);
                elseif cc(ii,mm)==2&&ii>mm
                    A(jj,kk,ii,mm) = A3(jj,kk);
                end
            end
        end
    end
end

A = shiftdim(A,2);


