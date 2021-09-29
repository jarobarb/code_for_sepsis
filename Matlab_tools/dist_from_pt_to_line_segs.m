function [mindists,dists,orthodists] = dist_from_pt_to_line_segs(ptsxs,ptsys,...
  linxs,linys)

  %  Points will correspond to row vectors, line segment points to column
  %  vectors
  npts = numel(ptsxs);
  nlin = numel(linxs);
  ptsxs = reshape(ptsxs,1,npts);
  ptsys = reshape(ptsys,1,npts);
  linxs = reshape(linxs,nlin,1);
  linys = reshape(linys,nlin,1);
  
  linxslower = linxs(1:end-1);
  linxsupper = linxs(2:end);
  linyslower = linys(1:end-1);
  linysupper = linys(2:end);

  tmpdists = sqrt(bsxfun(@minus,ptsxs,linxs).^2+...
    bsxfun(@minus,ptsys,linys).^2);
  dists = min(tmpdists(1:end-1,:),tmpdists(2:end,:));
  segdx = diff(linxs);
  segdy = diff(linys);
  lindist = sqrt(segdx.^2+segdy.^2);
  segdxnorm = segdx./lindist;
  segdynorm = segdy./lindist;
  pttoendpt1x = bsxfun(@minus,ptsxs,linxs(1:end-1));
  pttoendpt1y = bsxfun(@minus,ptsys,linys(1:end-1));
  dotdist = bsxfun(@times,segdxnorm,pttoendpt1x)+...
    bsxfun(@times,segdynorm,pttoendpt1y);
  orthodists = bsxfun(@times,segdxnorm,pttoendpt1y)-...
    bsxfun(@times,segdynorm,pttoendpt1x);
  inds = (dotdist >= 0) & bsxfun(@le,dotdist,lindist);
  dists(inds) = abs(orthodists(inds));
  mindists = min(dists);
  
%   mindists1 = mindists;
% 
%   dists = ptsxs';
%   orthodists = ptsxs';
%   for lc = 1:numel(linxs)-1
%     [dists(:,lc),orthodists(:,lc)] = dist_from_pt_to_line_seg(ptsxs,ptsys,...
%       [linxs(lc),linys(lc)],[linxs(lc+1),linys(lc+1)]);
%   end
%   dist = min(dists,[],2);
%   
%   keyboard

end