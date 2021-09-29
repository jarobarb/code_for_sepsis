function [v1M,v2M,fM] = read_fpde_table(varargin)

  filename = 'Table_exporting_and_importing_01.tbl';
    
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  fid = fopen(filename,'r');
  
  my_str = fgetl(fid);
  if ~isempty(strfind(my_str,'{'))
    while isempty(strfind(my_str,'}'))
      my_str = fgetl(fid);
    end
  end
  
  my_str = fscanf(fid,'%s',1);  
  nv1 = fscanf(fid,'%d',1);
  v1 = fscanf(fid,'%g',nv1);
  
  my_str = fscanf(fid,'%s',1);
  nv2 = fscanf(fid,'%d',1);
  v2 = fscanf(fid,'%g',nv2);
  
  [v1M,v2M] = ndgrid(v1,v2);
  
  my_str = fscanf(fid,'%s',2);
  my_str = fgetl(fid);
  my_vec = fscanf(fid,'%g');
  
  fM = reshape(my_vec,size(v1M));
  
  fclose(fid);

end