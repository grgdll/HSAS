function L2 = hsas_NIR_corr(L2, option)

# performs NIR correction for FICE 2022, 

if option == '1';           # (i) based on PML procesisng in Tilstone 2020

	L2.Rrs.data(:,201); # element 201 is at 750 nm
	Rrs_750_corr = ones(size(L2.Rrs.data, 1), length(L2.wv)).*L2.Rrs.data(:,201); 
	L2.Rrs.data = L2.Rrs.data - Rrs_750_corr; 

elseif option == '2';      # (ii) based on `simple NIR correction from hyperinspace - mimumum on [700,800] nm

	Rrs_minNIR_corr = ones(size(L2.Rrs.data, 1), length(L2.wv)).*min(L2.Rrs.data(:,176:226)')'; # 176 is 700 nm, 226 is 800 nm
	L2.Rrs.data = L2.Rrs.data - Rrs_minNIR_corr;
endif
