function [Natural_Frequencies,Modeshape]=modal_identification_fdd(AC,dt);
%Modal Identification 2017/06/22
%Author : Mohammad Nafisifard
%This Program using FDD Method for estimation of Modal properties. 
%Input:
%AC : Response data (rows = # of samples, columns = # number of channels)
%   dt :  Sampling period (sec)
%   REFs: 
%    [1] R. Brincker, L. Zhang, P. Andersen, "Modal Identification from 
%           Ambient Responses using Frequency Domain Decomposition", 18th 
%           International Modal Analysis Conference (IMAC XVIII), 2000. 
%    [2] S. Gade, N.B. Moller, H. Herlufsen, H. Konstantin-Hansen, 
%           "Frequency Domain Techniques for Operational Modal Analysis",
%           IMAC, Conference & Exposition on Structural Dynamics,
%           (IMAC-XXIV), 2006.

[n_row,n_col]=size(AC);
Num_points=2^(floor(log(length(AC))/log(2)));
AC_fft=fft(AC(1:Num_points,:));

% Singular Value Decomposition of the cross-spectrums
for kk=1:Num_points/2
    PSD_Ac=(conj(AC_fft(kk,:)))'*AC_fft(kk,:);
    [UU1,SS1,VV1]=svd((PSD_Ac(:,:)));
    ss1(kk)=abs(SS1(1,1));
    UU(:,kk)=UU1(:,1);
end

freq_domain=[0:(.5/dt)/(Num_points/2-1):.5/dt];

figure (100)
plot(freq_domain,log(ss1/2*length(freq_domain)));

%---------------------------------------------------------------------------------------------------
%... From plotted data, pick peaks on upper-most curve (dominant modes)

k=menu('# of PEAKS you identify in this plot?','1','2','3','4','5','6','7','8','9','10');
uiwait(msgbox('Click on PEAKS to extract Modal Parameters!','Left click on mouse','modal'));
[x1,y1] = ginput(k);
Natural_Frequencies=x1;

for jj=1:k
num_mode(jj)=round(x1(jj)/((.5/dt)/(Num_points/2-1))+1);
Modeshape_amp(:,jj)=UU(:,num_mode(jj));
phase_mode=phase(Modeshape_amp(:,jj));

for j1=1:n_col
    Phase_sign(j1,1)=sign(cos(phase_mode(1)-phase_mode(j1)));
end
    Modeshape(:,jj)=abs(Modeshape_amp(:,jj)).*Phase_sign;
end
