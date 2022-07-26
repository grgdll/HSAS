# HOW TO USE:  ls -d /Users/gdal/Cruises/AMT26/source/radiometry/hsas.source/output/L0/??? | xargs -P4  -n1 octave -q hsas_merge_L0.m
clear all

struct_levels_to_print (0);
#debug_on_error(1);

addpath ~/rad_functions

warning off

pkg load signal

VBS = false;  #VERBOSE

global main_dir
main_dir = "/Users/gdal/Cruises/AMT26/source/radiometry/hsas.source/";

fnin = argv;


doy2proc = "272";#strsplit(fnin{1}, "/"){end};


fn = glob([main_dir "output/L0/" doy2proc "/*L0.oct"]);

tic

# read first file to set the fields of L0
    tmp = load(fn{1});


# define size of L0 based on number of files and number of records in the first file

    L0.wv = tmp.L0.wv;

    n = length(tmp.L0.time)*4*length(fn);
    m = length(L0.wv);

    L0.time     = nan(n,1);
    L0.tilt     = nan(n,1);
    L0.roll     = nan(n,1);
    L0.pitch    = nan(n,1);
    L0.vaa_ths  = nan(n,1);

    L0.instr.Es.data = nan(n,m);
    L0.instr.Li.data = nan(n,m);
    L0.instr.Lt.data = nan(n,m);



# read all files and store in L0
irec = 1;
for ifn = 1:length(fn)

disp(fn{ifn});
fflush(stdout);

    tmp = load(fn{ifn});
    no = length(tmp.L0.time);

    L0.time(irec:irec+no-1)     = tmp.L0.time;
    L0.tilt(irec:irec+no-1)     = tmp.L0.tilt;
    L0.roll(irec:irec+no-1)     = tmp.L0.roll;
    L0.pitch(irec:irec+no-1)    = tmp.L0.pitch;
    L0.vaa_ths(irec:irec+no-1)  = tmp.L0.instr.Es.vaa_ths;

    L0.instr.Es.data(irec:irec+no-1,:) = tmp.L0.instr.Es.data;
    L0.instr.Li.data(irec:irec+no-1,:) = tmp.L0.instr.Li.data;
    L0.instr.Lt.data(irec:irec+no-1,:) = tmp.L0.instr.Lt.data;

    irec = irec+no;

endfor

# clear empty space in L0
    inan = find(isnan(L0.time));

    L0.time(inan)     = [];
    L0.tilt(inan)     = [];
    L0.roll(inan)     = [];
    L0.pitch(inan)    = [];
    L0.vaa_ths(inan)  = [];

    L0.instr.Es.data(inan,:) = [];
    L0.instr.Li.data(inan,:) = [];
    L0.instr.Lt.data(inan,:) = [];


toc

# save L0

    fnout = [main_dir "output/L0/hsas_" doy2proc "_L0.oct"];

    save("-binary", fnout, "L0");













