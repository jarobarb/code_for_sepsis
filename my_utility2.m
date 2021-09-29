function [vals,a,b] = my_utility2(varargin)

%   [a,b] = mapping('c1',0.5,'Bsource_in',1000,varargin{:});
  [a,b] = mapping_1d(varargin{:});
  fprintf('Healthy < Bsource = %g < Aseptic < Bsource = %g < Septic.\n',...
    a,b);
  
  %  All are normalized so that they are approximately between 0 and 1 when
  %  a and b are between 128 and 248
  %  Keeps the boundaries of the aseptic region near 128 and 248.  These
  %  two Bsources correspond to, we believe, septic and healthy outcomes in
  %  the actual experiments with the rats
  vals(1) = ((a-(128+1))^2+(b-(248-1))^2)/(2*((128+248)/2)^2);
  %  This keeps the aseptic region from going to zero
  vals(2) = (248-128)/(b-a);
  %  This keeps the aseptic region relatively centered
  vals(3) = ((b+a)/2-(128+248)/2)^2/(128-(128+248)/2)^2;
  
%   %  Changing c1, c2, c3 should make one of these more important than the
%   %  others.  We leave as the same importance for now
%   c1 = 1000; c2 = 1; c3 = 1;
%   val = c1*val1+c2*val2+c3*val3;
  
end