function my_struc = read_in_fpde_scalars(varargin)

  my_file_format = '*.p17';
  
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  ks = strfind(my_file_format,'/');
  if ~isempty(ks)
    my_root_dir = my_file_format(1:ks(end));
  else
    my_root_dir = './';
  end
  
  tmp_dir = dir(my_file_format);
  for j = 1:numel(tmp_dir)
    my_struc.name{j} = tmp_dir(j).name;
    tmp_struc = read_in_file([my_root_dir,my_struc.name{j}]);
    fns = fieldnames(tmp_struc);
    for fc = 1:numel(fns)
      my_struc.(fns{fc})(j) = tmp_struc.(fns{fc});
    end
  end

end

function my_struc = read_in_file(filename)

  fid = fopen(filename,'r');
  
  mystr = fgetl(fid);
  while ~strncmpi(mystr,'}',1)
    mystr = fgetl(fid);
  end
  while isempty(strfind(upper(mystr),'SUMMARY'))
    mystr = fgetl(fid);
  end
  
  mystr = fgetl(fid);
  while ischar(mystr)
    k = strfind(mystr,'=');
    tmp_name = mystr(1:k-1);
    for j = 1:2
      if isequal(tmp_name(1),' ')
        tmp_name = tmp_name(2:end);
      end
    end
    tmp_name = strrep(tmp_name,' ','_');
    tmp_name = strrep(tmp_name,'(','');
    tmp_name = strrep(tmp_name,')','');
    tmp_name = strrep(tmp_name,'*','');
    tmp_name = strrep(tmp_name,'-','m');
    my_struc.(tmp_name) = str2num(strrep(mystr(k+1:end),' ',''));
    mystr = fgetl(fid);
  end
  
  fclose(fid);

end