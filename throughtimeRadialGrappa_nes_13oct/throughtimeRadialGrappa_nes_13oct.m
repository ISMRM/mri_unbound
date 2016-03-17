function sigrec = throughtimeRadialGrappa_nes_13oct(data_u,acs,rseg,pseg,rep)
% 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
% To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US.
% By using this code, you implicitly agree to the terms in the license. 

% THROUGHTIME_RADIAL_GRAPPA reconstructs undersampled radial data using through-time
% radial GRAPPA, as described in Seiberlich, et al. MRM 2011 Feb;65(2):492-505.
%
% INPUTS
% data_u: undersampled k-space data [readout,projections,coils,measurements]
% acs: fully-sampled calibration data [readout,projections,coils,measurements]
% rseg: readout segment size
% pseg: projection segment size
% rep: number of measurement frames for calibration
% 
% OUTPUTS
% sigrec: reconstructed k-space [readout,projections,coils,measurements]
%
% EXAMPLE USAGE:
% sigrec = throughtime_radial_grappa(data_u,acs,8,4,26);
% -----------------------------------------------------------------------------------

%% ************* Error Checking ************** %%
assert(size(acs,4) >= rep, 'ACS does not have enough frames');
if rseg*pseg*rep <= 6*size(acs,3)
	warning('GRAPPA calibration is underdetermined. GRAPPA weights may not be stable.');
end

[nread,nproj,ncoils,~] = size(acs);
[~,nprojacc,~,nframes] = size(data_u);
accel = nproj/nprojacc; % acceleration factor


% At the borders of the auto-calibration data, the GRAPPA kernel will "overhang" k-space.
% To fix this, we pad the ACS data along the readout and projection dimensions by an amount
% "dr" and "dp". We have to use more padding when calibrating with a larger k-space segment
% (i.e. when rseg and pseg are larger).
RSEG = -floor(rseg/2) : floor( (rseg-1)/2 );
dr = max(abs(RSEG));
PSEG = -floor(pseg/2) : floor( (pseg-1)/2 );
dp = max(abs(PSEG));

dataref = zeros(nread+2*dr,nproj+2*dp+1,ncoils,rep,'single'); % dataref: the padded ACS data
dataref(dr+1:dr+nread,dp+1:dp+nproj,:,:) = acs;
dataref(dr+2:dr+nread,1:dp,:,:) = acs(end:-1:2,end-dp+1:end,:,:);
dataref(dr+2:dr+nread,dp+nproj+1:end,:,:) = acs(end:-1:2,1:dp+1,:,:);
dataref(1:dr,:,:,:) = dataref(2*dr:-1:dr+1,:,:,:);
dataref(dr+nread+1:end,:,:,:) = dataref(nread+1:dr+nread,:,:,:);
clear acs % clear ACS to save memory

% Initialize matrix to hold reconstructed k-space data. We repeat the first projection
% at the end so that the GRAPPA kernel never "overhangs" k-space.
sigrec = zeros(nread,nproj+1,ncoils,nframes,'single');
sigrec(:,1:accel:nproj,:,:) = data_u;
sigrec(end:-1:2,end,:,1:end-1) = sigrec(2:end,1,:,2:end);
sigrec(end:-1:2,end,:,end) = sigrec(2:end,1,:,end);
clear data_u % clear to save memory

kernelReps = rseg*pseg*rep;
overdetFctr = kernelReps/(6*ncoils);
fprintf('Through-time radial GRAPPA: Overdetermined factor = %.1f\n',overdetFctr)

% Initialize matrices for source and target points
src = zeros(6*ncoils,kernelReps);
targ = zeros((accel-1)*ncoils,kernelReps);

prevlen = 0;
for p = 1:accel:nproj
    msg = sprintf('Through-time radial GRAPPA: projection %d/%d\n',p,nproj);
    fprintf([repmat('\b',1,prevlen) '%s'],msg);
    prevlen = numel(msg);
    
    for r = 3:nread-2
        
        % calibration
        cnt = 1;
        for rr = 1:rseg
            for pp = 1:pseg
                rind = RSEG(rr)+r+dr;
                pind = PSEG(pp)+p+dp;
                src(:,cnt:cnt+rep-1) = reshape(dataref(rind-2:2:rind+2,pind:accel:pind+accel,:,:),6*ncoils,rep);
                targ(:,cnt:cnt+rep-1) = reshape(dataref(rind,pind+1:pind+accel-1,:,:),(accel-1)*ncoils,rep);
                cnt = cnt+rep;
            end
        end
        ws = pinv(src.')*targ.'; % GRAPPA weights for this particular location in k-space
        
        % reconstruction
        srcr = reshape(sigrec(r-2:2:r+2,p:accel:p+accel,:,:),6*ncoils,nframes); % source points from undersampled data
        targr = (srcr.'*ws).'; % target points
        sigrec(r,p+1:p+accel-1,:,:) = reshape(targr,[1,accel-1,ncoils,nframes]);
        
    end
end
sigrec = sigrec(:,1:nproj,:,:); % delete the extra projection at the end

msg = sprintf('Through-time radial GRAPPA: finished ;-)\n');
fprintf([repmat('\b',1,prevlen) '%s'],msg);