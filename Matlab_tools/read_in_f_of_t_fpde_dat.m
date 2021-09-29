function [ts,rest] = read_in_f_of_t_dat(varargin)

  filename = '34min_def/Temp_stats.dat';
  
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  fid = fopen(filename,'r');
  
  my_str = '';
  while isempty(strfind(my_str,'}'))
    my_str = fgetl(fid);
  end
  
  fscanf(fid,'%s',1);
  nt = fscanf(fid,'%d',1);
  
%   fprintf('%g\n',nt);
  
  ts = fscanf(fid,'%g',nt);
  
  while ~feof(fid)
    [~,b] = fscanf(fid,'%s',4);
    if b == 0, break; end
    fn = fscanf(fid,'%s',1);
    fscanf(fid,'%s',1);
    try
      rest.(fn) = fscanf(fid,'%g',nt);
    catch
      tmp_str = fgetl(fid);
      rest = fscanf(fid,'%g',nt);
    end
  end
  
  fclose(fid);

end