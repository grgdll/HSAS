function ths = hsas_rd_ths(fn)


#Instrument	InsFile	LogFile	Fitted	Immersed	Hex	CheckSum	TimeTags	F-TimeTags
#SATTHS0035	SATTHS0035b.tdf	2016-164-160729.raw	YES	NO	NO	NO	YES	YES

#Index	FRAME(COUNTER)	TIMER	COMP	PITCH	ROLL	DATETAG	TIMETAG2
#		sec	deg	deg	deg
#1	197	9598.0800000000	306.6000000000	-0.8000000000	-1.5000000000	2016-164	16:07:29.161
#2	198	9598.5600000000	306.8000000000	-0.6000000000	-1.7000000000	2016-164	16:07:29.640
#3	199	9599.0500000000	306.5000000000	-1.0000000000	-1.8000000000	2016-164	16:07:30.130
#4	200	9599.5300000000	306.4000000000	-1.1000000000	-1.7000000000	2016-164	16:07:30.698
#5	201	9600.1400000000	306.6000000000	-1.1000000000	-1.2000000000	2016-164	16:07:31.337

    DEBIAS = false;  # remove medfilt(angle, 101)

    ths.frame = [];
    ths.timer = [];
    ths.compass = [];
    ths.raw.pitch = [];
    ths.raw.roll = [];
    ths.time = [];


    for ifn = 1:length(fn)

        fid = fopen(fn{ifn});

            # skip header
                for irec = 1:5
                    fgets(fid);
                endfor

            fmt = "%f	%f	%f	%f	%f	%f	%f-%f	%f:%f:%f\n";
            tmp = fscanf(fid, fmt, [11, Inf])';

        fclose(fid);


        ths.frame = [ths.frame; tmp(:,2)];
        ths.timer = [ths.timer; tmp(:,3)];
        ths.compass = [ths.compass; tmp(:,4)];
        ths.raw.pitch = [ths.raw.pitch; tmp(:,5)];
        ths.raw.roll = [ths.raw.roll; tmp(:,6)];
        y0           = datenum(tmp(:,end-4),1,1)-1;
        ths.time     = [ths.time; y0+tmp(:,end-3)+tmp(:,end-2)/24+tmp(:,end-1)/24/60+tmp(:,end)/24/60/60];

        

    endfor


    
        ths.raw.tilt = hsas_cmp_tilt(ths.raw.roll, ths.raw.pitch);

    
    if DEBIAS

        w = 101;    

        ths.roll = ths.raw.roll - medfilt1(ths.raw.roll, w);
        ths.pitch = ths.raw.pitch - medfilt1(ths.raw.pitch, w);
        ths.tilt = ths.raw.tilt - medfilt1(ths.raw.tilt, w);

    else

        ths.roll = ths.raw.roll;
        ths.pitch = ths.raw.pitch;
        ths.tilt = ths.raw.tilt;

    endif





endfunction
