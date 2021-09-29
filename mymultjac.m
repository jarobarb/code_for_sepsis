function Js = mymultjac(t,ys,rps,ops,fps,strps)

  Js = zeros(numel(ys));
%   Js = sparse(numel(ys),numel(ys));

  for cc = 1:numel(rps)
    rp = rps(cc); op = ops(cc); fp = fps(cc); strp = strps(cc);
    inds = 4*(cc-1)+[1:5];

    Js(inds,inds) = myjac(t,ys(inds),rp,op,fp,strp);
  end

end