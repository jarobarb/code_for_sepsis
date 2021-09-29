function [xperave,f,inds,xwc,fftc,xperc,pdwrtpc] = ...
  my_pwelch(pres,NbrObs,nfft,noverlap,sF,varargin)

  %  Inputs
  %  pres-the data that we want to find the power spectral density for
  %  NbrObs-the size of the window to use in terms of number of data points
  %    when using overlapping windows for power spectral density analysis
  %  nfft-number of points to use for the discrete fourier transform
  %  noverlap-number of points to use for overlapping windows...0 no
  %    overlap...NbrObs-1 overlap all but 1 point (default is NbrObs/2)
  %  sF-sampling frequency...e.g. 5 Hz (5 cycles per second)
  %  different filtering windows can also be handed in (see my_window
  %    subfunction below)
  %
  %  Outputs
  %  xper-power associated with each frequency
  %  f-frequency associated with each power
  %  inds{nintervals}-indices used for each window
  %  xwc{nintervals}-corresponding windowed signal for a given window
  %  fftc{nintervals}-fft on just the window
  
  calc_pow_dens_wrt_period = false;
  window_type = 'Hann';
  
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  if isempty(nfft)
      nfft = max(256,2^nextpow2(NbrObs));
  end
  if isempty(noverlap)
      noverlap = round(NbrObs/2);
  end
  n = nargout;
% pres = detrend(pres);
% [Pxx,f] = pwelch(pres,NbrObs,[],[],sF);

  %  If last interval has length < NbrObs, include it (padded zeros added)
  nintervals = 1+ceil((numel(pres)-NbrObs)/(NbrObs-noverlap));
  %  If last interval has length < NbrObs, don't include it
%   s = warning('Using non-default behavior--ignoring incomplete samples',...
%     'backtrace','off');
%   disp(['Warning my_pwelch: ',...
%     'Using non-default behavior--ignoring incomplete samples.']);
  nintervals = 1+floor((numel(pres)-NbrObs)/(NbrObs-noverlap));
  intervals = round(bsxfun(@plus,[1,NbrObs],(NbrObs-noverlap)*[0:nintervals-1]'));
    
  hw = my_window(NbrObs,window_type);
  xper = zeros(nfft,nintervals);
  fftc = {};
  for k = 1:nintervals
    x1 = pres(intervals(k,1):min(intervals(k,2),numel(pres)));
    if length(x1) < NbrObs
      nintervals = nintervals-1;
      xper = xper(:,1:nintervals);
      intervals = intervals(1:nintervals,:);
      break;
    end
    x1(end+1:NbrObs) = 0;
    %  Detrending before applying window leads to undetrended windowed data
    %  which produces nonzero fft coefficients associated with period
    %  0...hence we start detrending after window (change 3/17/2015)
%     x1 = my_detrend(x1,'linear',[]);
    xw = x1.*hw;
    xw = my_detrend(x1.*hw,'linear',[]);
    tmp_fft = fft(xw,nfft);
    xper(:,k) = abs(tmp_fft).^2./norm(hw,2)^2;
    if n > 2
      inds{k} = intervals(k,1):min(intervals(k,2),numel(pres));
      xwc{k} = xw;
      fftc{k} = tmp_fft;
    else
      inds{1} = [];
      xwc{1} = [];
      fftc{1} = [];
    end
  end
  if n > 2
    xperc = xper./sF;
    xperc = xperc(1:size(xperc,1)/2+1,:);
    xperc(2:end-1,:) = 2*xperc(2:end-1,:);
  end
  xperave = mean(xper,2);
  xper = xper./sF;
  xperave = xperave./(sF);
  xper = xper(1:length(xperave)/2+1,:);
  xperave = xperave(1:length(xperave)/2+1);
  xper(2:end-1,:) = 2*xper(2:end-1,:);
  xperave(2:end-1) = 2*xperave(2:end-1);
  f = sF/nfft*(0:length(xperave)-1)';
  
  if calc_pow_dens_wrt_period
    per = 1./f;
    pow_dens_wrt_period = bsxfun(@rdivide,diff(cumtrapz(f,xper)),...
      abs(diff(per)));
    for k = 1:size(pow_dens_wrt_period,2)
      pdwrtpc{k} = pow_dens_wrt_period(:,k);
    end
  else
    pdwrtpc = [];
  end
%     max(abs(Pxx-xper))
end

function w = my_window(L,window_type)

  switch window_type
    case 'Hamming'
      w = 0.54-0.46*cos(2*pi*([0:L-1]')./(L-1));
    case 'Hann'
      w = (1/2)*(1-cos(2*pi*([0:L-1]')/(L-1)));
    otherwise
      w = ones(L,1);
  end

end