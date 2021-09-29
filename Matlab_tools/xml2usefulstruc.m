%  Given an structure from an xml file as produced by xml2struct and with
%  the xml file produced by PI Coresight, this changes the structure into a
%  more useful structure.
function my_struc = xml2usefulstruc(curr_struc)

  my_wrk_sheet = curr_struc.Workbook.Worksheet{2}.Table.Row;

  my_outer_str = '';
  j = 1;
  clear my_struc;
  while j < numel(my_wrk_sheet)
    if isfield(my_wrk_sheet{j}.Cell,'Data')
      %  Tells us this is an interesting row...i.e. bridge move times or
      %  underfeed/overfeed status level or resistance modifier level...
      %  Basically any variable that we had on our original SMART plot
      if ~isempty(strfind(my_wrk_sheet{j}.Cell.Data.Text,'-SMT-AF'))
        %  Get variable name
        my_var_name = my_wrk_sheet{j}.Cell.Data.Text(...
          strfind(my_wrk_sheet{j}.Cell.Data.Text,'|')+1:end);
        my_var_name = strrep(my_var_name,'.','p');
        %  Find start time
        j = j+1;
        day1 = datevec(strrep(my_wrk_sheet{j}.Cell{2}.Data.Text,'T',' '));
        %  Make sure we start at the beginning of the day
        day1(4:end) = 0;
        %  Skip end time
        j = j+1;
        %  Tell next loop to restart counting
        at_beginning = true;
      end
      j = j+1;
    else
      while ~isfield(my_wrk_sheet{j}.Cell,'Data')
        curr_datevec = datevec(strrep(my_wrk_sheet{j}.Cell{1}.Data.Text,...
          'T',' '));
        elapsed_time_in_sec = etime(curr_datevec,day1);
        if at_beginning
          myc = 1;
          at_beginning = false;
        else
          myc = myc+1;
        end
        my_struc.(my_var_name).t(myc) = elapsed_time_in_sec;
        my_struc.(my_var_name).val(myc) = str2num(...
          my_wrk_sheet{j}.Cell{2}.Data.Text);
        j = j+1;
      end
    end
  end

end
