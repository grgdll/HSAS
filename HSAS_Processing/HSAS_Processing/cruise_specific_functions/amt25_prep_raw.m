# divide AMT25 raw HSAS files into daily directories

din = "/Volumes/Public/DY110_Public/grg/AMT25/AMT25_raw/";

fnin = glob([din "*raw"]);


# find unique jday
doy = unique(cell2mat(fnin)(:,end-13:end-11),'rows');


for idoy = 1:length(doy)
	
	dtstr = datestr(datenum([2015 1 str2num(doy(idoy,:))]), "yyyymmdd")
	draw = [din "../Data/Underway_Rrs/" dtstr "/" dtstr "_raw"];
	mkdir(draw);
	
	system(["mv " din "2015-" doy(idoy,:) "*raw " draw]);
	
endfor









