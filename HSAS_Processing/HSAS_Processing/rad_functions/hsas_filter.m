function L2 = hsas_filter(L2, ikeep)
# set all L2 elements with indices different from ikeep to nan


n0 = length(L2.time);

f1 = fieldnames(L2);

for if1 = 1:length(f1)
   #disp(sprintf("f1=%s", f1{if1}))
   if strcmp(f1{if1},'instr')
      f2 = fieldnames( L2.(f1{if1}) );
      for if2 = 1:length(f2)# for each instrument
         #disp(sprintf("  f2=%s",f2{if2}))
         f3 = fieldnames(L2.(f1{if1}).(f2{if2}));
         for if3 = 1:length(f3)
            #disp(sprintf("  f3=%s",f3{if3}))
            if strcmp(f3{if3},"sn")
               disp(["Skipping " f3{if3}])
            else
               if length(L2.(f1{if1}).(f2{if2}).(f3{if3}))>length(ikeep) & ~isstr(L2.(f1{if1}).(f2{if2}).(f3{if3}))
                  L2.(f1{if1}).(f2{if2}).(f3{if3}) = L2.(f1{if1}).(f2{if2}).(f3{if3})(ikeep,:);
               else 
                  continue
               endif
            endif
         endfor
      endfor

   elseif strcmp(f1{if1},'gps') 
      f2 = fieldnames(L2.(f1{if1}));
      for if2 = 1:length(f2)
         L2.(f1{if1}).(f2{if2}) = L2.(f1{if1}).(f2{if2})(ikeep,:);
      endfor
      
   elseif strcmp(f1{if1},'oceanlogger') 
      f2 = fieldnames(L2.(f1{if1}));
      for if2 = 1:length(f2)
         L2.(f1{if1}).(f2{if2}) = L2.(f1{if1}).(f2{if2})(ikeep,:);
      endfor
      
   elseif strcmp(f1{if1},'hdr') | strcmp(f1{if1},'files')  
      L2.(f1{if1}) = L2.(f1{if1});
      
   elseif length(L2.(f1{if1}))==n0 & strcmp(f1{if1},'Lw')!=1 & strcmp(f1{if1},'doy')!=1
      L2.(f1{if1}) = L2.(f1{if1})(ikeep,:);
      
   elseif strcmp(f1{if1},'Lw') 
      L2.Lw.data = L2.Lw.data(ikeep,:);
      
   elseif strcmp(f1{if1},'Lwcorr') 
      L2.Lwcorr.data = L2.Lwcorr.data(ikeep,:);
      
#  else
#     continue
  endif
endfor

endfunction
