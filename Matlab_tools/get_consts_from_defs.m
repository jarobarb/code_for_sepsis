function get_consts_from_defs(varargin)

  filename = 'Stress_3d_Wedge_on_Side.pde';
  mat_file = '';
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  if ~isempty(mat_file)
    load(mat_file);
  end

  fid = fopen(filename,'r');

  my_str = fgetl(fid);

  if ~isempty(strfind(filename,'.pde'))
    while isempty(strfind(lower(my_str),'definitions'));
      my_str = fgetl(fid);
    end
  end

  while isempty(strfind(lower(my_str),'initial values')) &&...
      ~isequal(my_str,-1)
    my_str = strrep(my_str,'arctan','atan');
    my_str = strrep(my_str,'arccos','acos');
    my_str = strrep(my_str,'/','./');
    my_str = strrep(my_str,'*','.*');
    if isempty(strfind(my_str,'!'))% || ...
      %isempty(strfind(my_str,'{')) || ...
      %isempty(strfind(my_str,'}'))
      try, eval([my_str,';']); end
    else
      ind = strfind(my_str,'!');
      try, eval([my_str(1:ind-1),';']); end
    end
    my_str = fgetl(fid);
  end

  fclose(fid);
  
  clear varargin filename;
  save curr.mat
  
  display('Variables saved to curr.mat');

end