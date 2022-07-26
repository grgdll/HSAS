# divide AMT25 raw HSAS files into daily directories

din = "/Volumes/Public/DY110_Public/grg/AMT25/Data/Underway_Rrs/";

fnin = glob([din "2015????"]);


# find unique jday
dtstr = cell2mat(fnin)(:,end-7:end);


for idoy = 1:length(dtstr)
	
	
	dsat = [din dtstr(idoy,:) "/" dtstr(idoy,:) "_Satcon_extracted_raw"];
	mkdir(dsat);
	
	
	
endfor









