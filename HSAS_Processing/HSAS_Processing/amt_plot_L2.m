# plot L2 files to make sure everything is OK
clear all
close all

graphics_toolkit("gnuplot")

struct_levels_to_print (0);

warning off #Turn off warnings

addpath(strcat(pwd, "/cruise_specific_functions")) %FICE 2022: assumes code is run from ../HSAS_Processing 
addpath(strcat(pwd, "/rad_functions/"))
addpath(strcat(pwd, "/rad_functions/intwv"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman"))
addpath(strcat(pwd, "/rad_functions/DISTRIB_fQ_with_Raman/D_foQ_pa"))

# read input parameters for this cruise
input_parameters_hsas;


fnin = argv;
% fnin = {"20191017"};


disp('----------------------------------------------------------------------------------------------------------');
disp(fnin{1});
disp('----------------------------------------------------------------------------------------------------------');
fflush(stdout);


%DATESTR = fnin{1};
DATESTR = strsplit(fnin{1},'/'){end-1};

XLIMS = [400 750];

DIN_L2 = strrep(DIN_L1,'L1','L2')
fn = glob([DIN_L2 DATESTR "/*mat"]);

DOUT = [DIN_L2 DATESTR "/plot/"];

if ~isdir(DOUT)
	mkdir(DOUT);
endif



for ifn =1:length(fn)

   clear L2 L1_f

   load(fn{ifn});

   # select range of wv to be plotted
   iwv = find(L2.wv>=XLIMS(1) & L2.wv<=XLIMS(2));

   figure(1, 'visible', 'off')

   clf
      subplot(241)
         hold on
            plot(L1_f.time, L1_f.instr.Es.data(:,L1_f.wv==550), "ko-")
            plot(L2.time, L2.instr.Es.data(:,L2.wv==550), "ro")
         datetick("x", "HH");
         xlabel ("hour")
         ylabel ("Es(550) [uW/cm^2/nm]")

      subplot(242)
         hold on
            plot(L1_f.time, L1_f.phi, "ko-")
            plot(L2.time, L2.phi, "ro")
         datetick("x", "HH");
         xlabel ("hour")
         ylabel ("\Delta Azimuth(sensor-sun) [degrees]")

      subplot(243)
         hold on
            plot(L2.wv(iwv), L2.Rrs.data(:,iwv), "r-")
            plot(L2.wv(iwv), L2.wv(iwv)*0, 'k')
            plot(L2.wv(iwv), L2.Rrs.data(L2.phi<90,iwv), 'k')
         xlabel("wavelength [nm]")
         ylabel("Rrs [1/sr]")
         set(gca, 'pos', [0.55699   0.58682   0.30   0.33818])

         ylim([-0.0025 0.020])


      #subplot(244)


      subplot(245)
         hold on
            plot(L2.wv(iwv), L2.instr.Lt.data(:,iwv), "r-")
         xlabel("wavelength [nm]")
         ylabel("Lt [uW/cm^2/nm/sr]")

      subplot(246)
         hold on  
#            plot(L1_f.wv(iwv), L1_f.instr.Lt.data, "k-")
            plot(L2.wv(iwv), L2.instr.Es.data(:,iwv), "r-")
         xlabel("wavelength [nm]")
         ylabel("Es [uW/cm^2/nm]")

      subplot(247)
         hold on
#            plot(L1_f.wv(iwv), L1_f.instr.Li.data, "k-")
            plot(L2.wv(iwv), L2.instr.Li.data(:,iwv), "r-")
         xlabel("wavelength [nm]")
         ylabel("Li [uW/cm^2/nm/sr]")

      #subplot(248)
      #   hold on
      #      plot(L1_f.time, L1_f.instr.Lt.data(:,L1_f.wv==440), "ko-")
      #      plot(L2.time, L2.instr.Lt.data(:,L2.wv==440), "ro")
      #   datetick("x", "HH")
      #   xlabel("hour")
      #   ylabel("Lt(440) [uW/cm^2/nm/sr]")
      subplot(248)
         hold on
            #plot(L1_f.time, L1_f.instr.Lt.data(:,L1_f.wv==440), "ko-")
         plot(L2.time, pi*L2.instr.Li.data(:,L2.wv==400)./L2.instr.Es.data(:,L2.wv==400), "bo;400 nm;")
         plot(L2.time, pi*L2.instr.Li.data(:,L2.wv==650)./L2.instr.Es.data(:,L2.wv==650), "ro;650 nm;")
         datetick("x", "HH")
         xlabel("hour")
         ylabel("\pi*Li(\lambda)/Es(\lambda) [-]")

set(gcf, 'paperposition', [0.25 0.25 14 6])

fnout = strrep(fn{ifn}, "mat", "png");

fnout = strrep(fnout, DATESTR, [DATESTR "/plot"]);
print(fnout, "-dpng")

endfor



