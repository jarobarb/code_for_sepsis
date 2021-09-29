function [dist,orthodist] = dist_from_pt_to_line_seg(ptxs,ptys,...
  linendpt1,linendpt2)

  dist = min(sqrt((ptxs-linendpt1(1)).^2+(ptys-linendpt1(2)).^2),...
    sqrt((ptxs-linendpt2(1)).^2+(ptys-linendpt2(2)).^2));
  linvec = linendpt2-linendpt1;
  lindist = norm(linvec);
  linnorm = linvec./lindist;
  pttoendpt1x = ptxs-linendpt1(1);
  pttoendpt1y = ptys-linendpt1(2);
  dotdist = linnorm(1)*pttoendpt1x+linnorm(2)*pttoendpt1y;
  orthodist = (linnorm(1)*pttoendpt1y-linnorm(2)*pttoendpt1x);
  inds = (dotdist >= 0) & (dotdist <= lindist);
  dist(inds) = abs(orthodist(inds));

end