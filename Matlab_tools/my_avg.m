function avg_vals = my_avg(vals,n)

  avg_vals = vals(:);
  for ac = 1:n
    avg_vals = [avg_vals(1);...
      (avg_vals(1:end-2)+2*avg_vals(2:end-1)+avg_vals(3:end))./4;...
      avg_vals(end)];
  end
  avg_vals = reshape(avg_vals,size(vals));

end