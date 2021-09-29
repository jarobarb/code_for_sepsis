function convert_xy_to_logs(hl,hp,lower_bound)
  
  for hh = [hl;hp]'
    xs = max(get(hh,'XData'),lower_bound);
    set(hh,'XData',log10(xs));
    ys = max(get(hh,'YData'),lower_bound);
    set(hh,'YData',log10(ys));
  end
%   keyboard

end