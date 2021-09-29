function write_fpde_table(x,y,f,varargin)

  out_file = 'temp.tbl';
  
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  xs = unique(x);
  ys = unique(y);
  
  fid = fopen(out_file,'w');
  fprintf(fid,'\nX %d\n',numel(xs));
  fprintf(fid,'%12.10g %12.10g %12.10g %12.10g %12.10g\n',xs);
  fprintf(fid,'\nY %d\n',numel(ys));
  fprintf(fid,'%12.10g %12.10g %12.10g %12.10g %12.10g\n',ys);
  
  fprintf(fid,'\nDATA {data}\n');
  for yc = 1:numel(ys)
    fprintf(fid,'%12.10g %12.10g %12.10g %12.10g %12.10g\n',f(:,yc));
    fprintf(fid,'\n');
  end
  
  fclose(fid);
  
end