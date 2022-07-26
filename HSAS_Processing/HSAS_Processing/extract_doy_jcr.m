# this script separates the seatex-htg.ACO file into files containing different doy 
ship_data_stream = argv(){1};

DIN = "/data/datasets/cruise_data/active/AMT24/Ship_data/Compress/Compress/";
fn = ["/data/datasets/cruise_data/active/AMT24/Ship_data/Compress/Compress/" ship_data_stream ".ACO"];

DOUT = [DIN 'days/' ship_data_stream '/'];   
if ~exist(DOUT)
   mkdir(DOUT);
endif

#id = load(fn);
for idoy =265:306
	disp(num2str(idoy));
	system(	sprintf(['grep 2014,%3u   %s > %sdays/' ship_data_stream '/' ship_data_stream '.%3u' ], idoy, fn, DIN, idoy)  );
	fflush(stdout);	
 	
endfor


