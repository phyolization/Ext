function [SW, SigLam] = steeringWeight_Loss(sigma,varargin)
%STEERINGWEIGHT calculates the steering weight of an assemblage
% This function has one required argument:
%  sigma: a 4-D array, containing the members of the assemblage. The first 
%   two dimensions contain the (unnormalised) quantum states, while the
%   remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%   \sigma_a|x. 
%
% SW = steeringWeight(sigma) returns the steering weight SW of the
% assemblage sigma.
%
% This function has one optional argument:
%   He: (default 1)
%
% ->requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
% ->authors: Paul Skrzypczyk, Daniel Cavalcanti (Modified by Qiang Z)
% ->last updated: Jan 10, 2021

[He,of] = opt_args({1,1},varargin{:});
% if unspecified, it is assumed the standard steering weight is required.
% He: Heralding efficiency
% of: Optimization factor

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice
pri=oa+1; %the priori scenario
Ndet = (pri)^ma; %number of deterministic behaviours

SingleParty = zeros(pri,ma,Ndet); % initialise array
for lam = 0:Ndet-1
    lamdec = dec2base(lam,pri,ma)-'0'; % generates the string of outcomes a 
                                %(for each x), for the given variable lam
                                %将每个字符的ascii码减去‘0’的ascii码即获得该字符串的数列
    for x = 1:ma
        for a = 0:pri-1
;           SingleParty(a+1,x,lam+1) = (lamdec(x) == a);...   +(lamdec(x) == a).*(1-He); %
            % % probability = 1 if a = lamdec(x), 0 otherwise
        end
    end
end
% generate array containing the single party distributions

SingleParty_N = zeros(dB,dB,Ndet,pri*ma);
for i = 1:ma
    for j = 1:pri
        for k = 1:pri^ma
          SingleParty_N(1:dB,1:dB,k,j,i) = SingleParty(j,i,k);
        end
    end
end


% check that the assemblage is valid
if  NSAssemblage(sigma) == 0
    error('assemblage is not valid')
end

sigR = sum(sigma(:,:,:,1),3);
Sigma=sigma.*He;
Sigma(:,:,pri,:)=repmat(sigR.*(1-He),[1,1,1,ma]);
% sigR = reduced state of Bob, necessary for consistent steering robust.

% NOTE: Here we use the dual formulation of the steering weight.

cvx_begin sdp quiet
    
    variable SigLam(dB,dB,Ndet) hermitian semidefinite
    % members of the steering functional

    maximise real(sum(reshape(SigLam.*repmat(eye(dB),[1,1,Ndet]),1,[])))
    % sum(reshape(X.*conj(Y),1)) equals to the Trace!!! (if Y is hermitian)
    % Fascinating!!!
    
    subject to
 
    for i = 1:ma
        for j=1:oa

        - sum(SingleParty_N(:,:,:,j,i).*SigLam,3) ...
             + 1.000001*Sigma(:,:,j,i) == hermitian_semidefinite(dB); %
                  
        end

        - sum(SingleParty_N(:,:,:,pri,i).*SigLam,3) ...
             + 1.000001*Sigma(:,:,pri,i) == hermitian_semidefinite(dB); %
     
    end
    
cvx_end

SW = 1-cvx_optval;
    
end
