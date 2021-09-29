function [ccell,cxs,cys] = extract_contours(cmatrix)

  start = 1; myc = 1;
  while start < size(cmatrix,2)
    npts = cmatrix(2,start);
    ccell{myc} = cmatrix(:,start+[1:npts]);
    cxs{myc} = ccell{myc}(1,:);
    cys{myc} = ccell{myc}(2,:);
    start = start+npts+1; myc = myc+1;
  end

end