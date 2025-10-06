%FEVD - Compute the forecast error variance decomposition
%
%
% uses [hx^0, hx^1, ....]eta*eta'*[hx^0,hx^1,...]'


function [outY,outX, outT] = fevd(gx,hx,eta,h)

%Dimensions
nh = length(h);
nx = size(hx,1);
ny = size(gx,1);
neps = size(eta,2);

%[hx^0, hx^1, ....]
H = zeros(nx,nh*nx);
H(:,1:nx) = eye(nx);
for hh = 2:max(h)
   H(:,(hh-1)*nx+(1:nx)) = hx^(hh-1); 
end


outY = zeros(neps,ny,nh);
outX = zeros(neps,nx,nh);
outT = zeros(1,ny+nx,nh);
hcount = 1;
for hh = h
    
    %Total FEV
    ETA = repmat_diag(eta,hh);
    sigX = H(:,1:hh*nx)*(ETA*ETA')*H(:,1:hh*nx)';
    sigY = gx*sigX*gx';
    outT(:,:,hcount) = [diag(sigY)',diag(sigX)'];
    %Taking shocks one at a time
    for jj = 1:neps
        ETA = repmat_diag(eta(:,jj),hh);
        sigX_part = H(:,1:hh*nx)*(ETA*ETA')*H(:,1:hh*nx)';
        sigY_part = gx*sigX_part*gx';
        outY(jj,:,hcount) = diag(sigY_part)'./diag(sigY)';
        outX(jj,:,hcount) = diag(sigX_part)'./diag(sigX)';
    end
    
    hcount = hcount+1;
    
end