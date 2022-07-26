% Function TRUEW.M  
% Converted from C true wind computation code truewind.c .

% Created: 12/17/96
% Comments updated: 10/01/97
% Developed by: Shawn R. Smith and Mark A. Bourassa
% Programmed by: Mylene Remigio
% Last updated:  9/30/2014
% Direct questions to:  wocemet@coaps.fsu.edu 
% 
% 9/30/2014 : If the true wind has speed and its coming from the north 
%	     then its direction should be 360deg. The problem was fixed
%	     in which the programs’s output showed 0deg instead of 360deg.
%
% This routine will compute meteorological true winds (direction from
% which wind is blowing, relative to true north; and speed relative to
% the fixed earth).
%
% INPUT VALUES:
%
% num	integer		Number of observations in input (crse, cspd,
%                       wdir, wspd, hd) and output (adir, tdir, tspd)
%		        data arrays. ALL ARRAYS MUST BE OF EQUAL LENGTH.
% sel	integer		Sets option for diagnostic output.  There are
%			four settings:
%
%			Option 4:  Calculates true winds from input
%				   arrays with no diagnostic output or
%				   warnings. NOT RECOMMENDED.
%			Option 3:  [DEFAULT] Diagnostic output lists the
%                                  array index and corresponding 
%				   variables that either violate the 
%    				   range checks or are equal to the 
%				   missing value. An additional table
%                                  lists the number of observation times
%                                  with no missing values, some (but not
%				   all)  missing values, and all missing
%				   values; as well as similar totals for
%				   the observation times that fail the 
%				   range checks. Range checks identify 
%				   negative input values and verify 
%				   directions to be between 0 and 360 
%				   degrees.
%			Option 2:  In addition to the default 
%				   diagnostics (option 3), a table of 
%				   all input and output values for
%				   observation times with missing data
%                                  is provided.
%			Option 1:  Full diagnostics -- In addition to
%				   the diagnostics provided by option 2 
%				   and 3, a complete data chart is
%                                  output. The table contains input and
%                                  output values for all observation 
%                                  times passed to truewind.
%
% crse	real array	Course TOWARD WHICH the vessel is moving over
%			the ground. Referenced to true north and the 
%                       fixed earth.
% cspd	real array	Speed of vessel over the ground. Referenced
%			to the fixed earth.
% hd	real array	Heading toward which bow of vessel is pointing.
%			Referenced to true north.
% zlr	real		Zero line reference -- angle between bow and
%			zero line on anemometer.  Direction is clockwise
%			from the bow.  (Use bow=0 degrees as default 
%			when reference not known.)
% wdir 	real array	Wind direction measured by anemometer,
%			referenced to the ship.
% wspd	real array	Wind speed measured by anemometer,referenced to
%			the vessel's frame of reference.
% wmis	real array	Five element array containing missing values for
%			crse, cspd, wdir, wspd, and hd. In the output, 
%                       the missing value for tdir is identical to the 
%                       missing value specified in wmis for wdir. 
%                       Similarly, tspd uses the missing value assigned 
%			to wmis for wspd.

% *** WDIR MUST BE METEOROLOGICAL (DIRECTION FROM)! CRSE AND CSPD MUST 
%     BE RELATIVE TO A FIXED EARTH! ***

% OUTPUT VALUES:

% tdir	real array	True wind direction - referenced to true north
%                       and the fixed earth with a direction from which 
%			the wind is blowing (meteorological).
% tspd	real array	True wind speed - referenced to the fixed earth.
% adir	real array	Apparent wind direction (direction measured by
%			wind vane, relative to true north). IS 
%                       REFERENCED TO TRUE NORTH & IS DIRECTION FROM
%                       WHICH THE WIND IS BLOWING. Apparent wind
%			direction is the sum of the ship relative wind
%                       direction (measured by wind vane relative to the
%                       bow), the ship's heading, and the zero-line
%			reference angle.  NOTE:  The apparent wind speed
%			has a magnitude equal to the wind speed measured
%		  	by the anemometer.

% DIAGNOSTIC OUTPUT:

% nw	integer		Number of observation times for which tdir and
%                       tspd were calculated (without missing values)
% nwpm	integer		Number of observation times with some values
%  			(crse, cspd, wdir, wspd, hd) missing.  tdir, 
%			tspd set to missing value. 
% nwam	integer		Number of observation times with all values
%  			(crse, cspd, wdir, wspd, hd) missing. tdir, 
%			tspd set to missing value.
% nwf	integer		Number of observation times where the program
%			fails -- at least one of the values (crse, cspd,
%			wdir, wspd, hd) is invalid 

%************************************************************************

function truew(num, sel, crse, cspd, wdir, zlr, hd, adir, wspd, ...
    wmis, tdir, tspd, nw, nwam, nwpm, nwf)

% INITIALIZE VARIABLES

x = 0;
y = 0;
nw = 0;
nwam = 0;
nwpm = 0;
nwf = 0;
dtor = pi/180;		%degrees to radians conversion

% LOOP OVER 'NUM' VALUES
i = 1;
while( i <= num)

%   CHECK COURSE, SHIP SPEED, HEADING, WIND DIRECTION, AND WIND 
%   SPEED FOR VALID VALUES (i.e. neither missing nor outside 
%   physically acceptable ranges).

    if( ( ((crse(i) < 0) || (crse(i)> 360)) && (crse(i) ~= wmis(1)) ) || ...
        ( (cspd(i) < 0) && (cspd(i) ~= wmis(2)) )|| ...
        ( ((wdir(i) < 0) || (wdir(i) > 360)) && (wdir(i) ~= wmis(3)) ) || ...
        ( (wspd(i) < 0) && (wspd(i) ~= wmis(4)) ) || ...
        ( ((hd(i) < 0) || (hd(i) > 360)) && (hd(i) ~= wmis(5)) ) ) 
%	WHEN SOME OR ALL OF INPUT DATA FAILS RANGE CHECK, TRUE WINDS ARE SET
%	TO MISSING. STEP INDEX FOR INPUT VALUE(S) BEING OUT OF RANGE   
	    nwf = nwf + 1;
        tdir(i) = wmis(3);
	    tspd(i) = wmis(4);   


        if( (crse(i) ~= wmis(1)) && (cspd(i) ~= wmis(2)) && ... 
	  (wdir(i) ~= wmis(3)) && (wspd(i) ~= wmis(4)) && ...
	  (hd(i) ~= wmis(5)) )
%       STEP INDEX FOR ALL INPUT VALUES BEING NON-MISSING  
	  nw = nw + 1;
        else 
                if( (crse(i) ~= wmis(1)) || (cspd(i) ~= wmis(2)) || ...
          	(wdir(i) ~= wmis(3)) || (wspd(i) ~= wmis(4)) || ...
          	(hd(i) ~= wmis(5)) )
%       STEP INDEX FOR PART OF INPUT VALUES BEING MISSING 
                   nwpm = nwpm + 1;
                else
%       STEP INDEX FOR ALL INPUT VALUES BEING MISSING 
                   nwam = nwam + 1;
                end
        end

%   WHEN COURSE, SHIP SPEED, HEADING, WIND DIRECTION, AND WIND 
%   SPEED ARE ALL IN RANGE AND NON-MISSING, THEN COMPUTE TRUE WINDS.

    elseif( (crse(i) ~= wmis(1)) && (cspd(i) ~= wmis(2)) && ...
            (wdir(i) ~= wmis(3)) && (wspd(i) ~= wmis(4)) && ...
	  (hd(i) ~= wmis(5)) )
    	nw = nw + 1;

%	CONVERT FROM NAVIGATIONAL COORDINATES TO ANGLES COMMONLY USED 
% 	IN MATHEMATICS
    	mcrse = 90 - crse(i);

%         KEEP THE VALUE BETWEEN 0 AND 360 DEGREES  
	    if( mcrse <= 0.0 ) 
		    mcrse = mcrse + 360.0;
	    end

%	CHECK ZLR FOR VALID VALUE.  IF NOT VALID, SET EQUAL TO ZERO.
	    if( (zlr < 0.0) || (zlr > 360.0) ) 
		    zlr = 0.0;
	    end

%	CALCULATE APPARENT WIND DIRECTION 
	     adir(i) = hd(i) + wdir(i) + zlr;

%	KEEP ADIR BETWEEN 0 AND 360 DEGREES 
	    while( adir(i) >= 360.0 ) 
		    adir(i) = adir(i) - 360.0;	 
	    end

%         CONVERT FROM METEOROLOGICAL COORDINATES TO ANGLES COMMONLY USED   
%         IN MATHEMATICS 
	    mwdir = 270.0 - adir(i);

%         KEEP MDIR BETWEEN 0 AND 360 DEGREES 
          if( mwdir <= 0.0 ) 
	    mwdir = mwdir + 360.0;
          end
          if( mwdir > 360.0 )
	    mwdir = mwdir - 360.0;
          end

%        DETERMINED THE EAST-WEST VECTOR COMPONENT AND THE NORTH-SOUTH
%        VECTOR COMPONENT OF THE TRUE WIND 
	x = (wspd(i) .* cos(mwdir .* dtor)) + (cspd(i)  .* cos(mcrse .* dtor));
	y = ((wspd(i) .* sin(mwdir .* dtor)) + (cspd(i) .* sin(mcrse .* dtor)));

%        USE THE TWO VECTOR COMPONENTS TO CALCULATE THE TRUE WIND SPEED  
	tspd(i) = sqrt((x .* x) + (y .* y));
	calm_flag = 1;

%        DETERMINE THE ANGLE FOR THE TRUE WIND 
	if(abs(x) > 0.00001) 
	    mtdir = (atan2(y,x)) / dtor;
          else
	    if(abs(y) > 0.00001) 
	        mtdir = 180.0 - (90.0 * y) / abs(y);
	    else
%               THE TRUE WIND SPEED IS ESSENTIALLY ZERO: WINDS ARE CALM
%               AND DIRECTION IS NOT WELL DEFINED   
                  mtdir = 270.0;
                  calm_flag = 0;
	    end
	end

%        CONVERT FROM THE COMMON MATHEMATICAL ANGLE COORDINATE TO THE 
%        METEOROLOGICAL WIND DIRECTION  
	tdir(i) = 270.0 - mtdir;

%        MAKE SURE THAT THE TRUE WIND ANGLE IS BETWEEN 0 AND 360 DEGREES 
	while(tdir(i) < 0.0) 
	    tdir(i) = (tdir(i) + 360.0) .* calm_flag;
	end
	while(tdir(i) > 360.0) 
	    tdir(i) = (tdir(i) - 360.0) .* calm_flag;
	end
%        ENSURE WMO CONVENTION FOR TDIR=360 FOR WIN FROM NORTH AND TSPD > 0 
	if (calm_flag == 1 && (tdir(i) < 0.0001))
              tdir(i) = 360.0;
	end
	x = 0.0; 
	y = 0.0;
    



    else
    
    
        if( (crse(i) ~= wmis(1)) || (cspd(i) ~= wmis(2)) || ...
            (wdir(i) ~= wmis(3)) || (wspd(i) ~= wmis(4)) || ...
            (hd(i) ~= wmis(5)) )
	nwpm = nwpm + 1;
          tdir(i) = wmis(3);
          tspd(i) = wmis(4);
	
%         WHEN COURSE, SHIP SPEED, APPARENT DIRECTION, AND WIND SPEED 
%         ARE ALL IN RANGE BUT ALL OF THESE INPUT VALUES ARE MISSING,
%         THEN SET TRUE WIND DIRECTION AND SPEED TO MISSING.      

	else
              nwam = nwam + 1;
              tdir(i) = wmis(3);
              tspd(i) = wmis(4);
	end


    end

%   INCREASE COUNTER
    i = i + 1;
    end

%   OUTPUT SELCTION PROCESS
    switch(sel)
        case 1
            full(num, crse, cspd, wdir, zlr, hd, adir, wspd, tdir, tspd), ...
            missing_values(num, crse, cspd, wdir, hd, wspd, tdir, tspd, wmis), ...
            truerr(num, crse, cspd, hd, wdir, wspd, wmis, nw, nwpm, nwam, nwf);
        case 2
            missing_values(num, crse, cspd, wdir, hd, wspd, tdir, tspd, wmis), ...
            truerr(num, crse, cspd, hd, wdir, wspd, wmis, nw, nwpm, nwam, nwf);
        case 3
	  truerr(num, crse, cspd, hd, wdir, wspd, wmis, nw, nwpm, nwam, nwf);
        otherwise
	  fprintf('Selection not valid. Using selection #3 by default. \n'), ...
	  truerr(num, crse, cspd, hd, wdir, wspd, wmis, nw, nwpm, nwam, nwf);
end


end


% **********************************************************************
%                          OUTPUT SUBROUTINES 
% **********************************************************************

%     Function:  FULL
%      Purpose:  Produces a complete data table with all values. 
%                Accessed only when selection #1 is chosen.


function full(num, crse, cspd, wdir, zlr, hd, adir, wspd, tdir, tspd)
fprintf('------------------------------------------------------------------------------------\n');
fprintf('\n                                   FULL TABLE\n');
fprintf('                                  ************\n');
fprintf('  index  course  sspeed  windir  zeroln  shiphd |  appspd |  appdir  trudir  truspd\n')
for j = 1:num

      fprintf('%7d %7.1f %7.1f %7.1f %7.1f %7.1f | %7.1f | %7.1f %7.1f %7.1f\n', ...
          (j),crse(j),cspd(j),wdir(j),zlr,hd(j),wspd(j),adir(j),tdir(j),tspd(j));
end


fprintf('\n                   NOTE:  Wind speed measured by anemometer is identical\n');
fprintf('                          to apparent wind speed (appspd).\n');
fprintf('\n------------------------------------------------------------------------------------\n');
end

% **********************************************************************

%    Function:  MISSING_VALUES
%     Purpose:  Produces a data table of the data with missing values. 
%               Accessed when selection #1 or #2 is chosen.


function missing_values(num, crse, cspd, wdir, hd, wspd, tdir, tspd, wmis)
fprintf('\n                               MISSING DATA TABLE\n');
fprintf('                              ********************\n');
fprintf('          index  course  sspeed  windir  shiphd  appspd  trudir  truspd\n');

for j = 1:num

      if( (crse(j) ~= wmis(1)) && (cspd(j) ~= wmis(2)) &&  ...
          (wdir(j) ~= wmis(3)) && (wspd(j) ~= wmis(4)) &&  ...
          (hd(j) ~= wmis(5)) )
        continue;
      else                   
         fprintf('        %7d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n', ...
                (j), crse(j), cspd(j), wdir(j), hd(j), wspd(j), tdir(j), tspd(j));
      end
end

   fprintf('\n------------------------------------------------------------------------------------\n');
end   

% **********************************************************************

%    Function:  TRUERR
%     Purpose:  List of where range tests fail and where values are
%               invalid.  Also prints out number of records which are
%               complete, incomplete partially, incomplete entirely, and
%               where range tests fail.  Accessed when selection #1, #2,
%               #3, or the default is chosen.


 
function truerr(num, crse, cspd, hd, wdir, wspd, wmis, nw, nwpm, nwam, nwf)
   fprintf('\n                               TRUEWINDS ERRORS\n');
   fprintf('                              ******************\n');


   for i=1:num                    
      if( ((crse(i) < 0) || (crse(i) > 360)) && (crse(i) ~= wmis(1)) )                                             
        fprintf('        Truewinds range test failed.  Course value #%d invalid.\n',(i));
      end
      if( (cspd(i) < 0) && (cspd(i) ~= wmis(2)) )
        fprintf('        Truewinds range test failed.  Vessel speed value #%d invalid.\n',(i));
      end
      if( ((wdir(i) < 0) || (wdir(i) > 360)) && (wdir(i) ~= wmis(3)) )     
        fprintf('        Truewinds range test failed.  Wind direction value #%d invalid.\n',(i));
      end
      if( (wspd(i) < 0) && (wspd(i) ~= wmis(4)) ) 
        fprintf('        Truewinds range test failed.  Wind speed value #%d invalid.\n',(i));
      end
      if( ((hd(i) < 0) || (hd(i) > 360)) && (hd(i) ~= wmis(5)) )                                             
        fprintf('        Truewinds range test failed.  Ship heading value #%d invalid.\n',(i));
      end
   end

   fprintf('\n');

   for i =1:num 
      if(crse(i) == wmis(1))      
        fprintf('        Truewinds data test:  Course value #%d missing.\n', (i));
      end
      if(cspd(i) == wmis(2))
        fprintf('        Truewinds data test:  Vessel speed value #%d missing.\n', (i));
      end
      if(wdir(i) == wmis(3))
        fprintf('        Truewinds data test:  Wind direction value #%d missing.\n', (i));
      end
      if(wspd(i) == wmis(4))
        fprintf('        Truewinds data test:  Wind speed value #%d missing.\n', (i));
      end
      if(hd(i) == wmis(5))
        fprintf('        Truewinds data test:  Ship heading value #%d missing.\n', (i));
      end
   end    


   fprintf('\n------------------------------------------------------------------------------------\n');
   fprintf('\n                                 DATA REVIEW\n');
   fprintf('                                *************\n');
   fprintf('                            no data missing = %4d\n', nw);
   fprintf('                       part of data missing = %4d\n', nwpm);
   fprintf('                           all data missing = %4d\n', nwam);
   fprintf('                         failed range tests = %4d\n', nwf);    
end
   
   
     
   
   
   
   
   
   
