function [k,klabels] = define_default_ks

  %  In an effort to not have multiple definitions in multiple locations
  %  (making it so that if we decide to change parameter values, we have
  %  change it in all of the locations where those parameter definitions
  %  are found), we restored this program to its original intent.  We note
  %  that many programs still overwrite these values, especially when
  %  exploring parameter space.

  %  See my_sensitivity_script7.m, section "Alternate optimization for "kA"
  %  and "k5" set.
  k = [5.0046602934939192764
         70.999380840544219495
       0.030009711419666851295
         80.838762887398843304
         1468.9689068779709942
          7.1834182161733215466]';
  
  klabels = {'B_infy','f','kD1','Bmf','nu4','kA'};
  
end