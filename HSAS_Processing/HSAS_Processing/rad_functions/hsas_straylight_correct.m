function data = hsas_straylight_correct(data)
input_parameters_hsas;

VBS = false;

#----Read Straylight Distribution Matrix
#----_Es:
fn=[DIR_SLCORR,FN_SLCORR_ES];
D_Es = load(fn);
D_Es = D_Es / norm(D_Es);

#----_Li:
fn=[DIR_SLCORR,FN_SLCORR_LI];
D_Li = load(fn);
D_Li = D_Li / norm(D_Li);

#----_Lt:
fn=[DIR_SLCORR,FN_SLCORR_LT];
D_Lt = load(fn);
D_Lt = D_Lt / norm(D_Lt);

wv = data.wv;
I = eye(length(wv),length(wv));

#----Correction for Es
A = I +D_Es;
disp('** Straylight correction for Es,wait!')
for id = 1:length(data.instr.Es.data(:,1));
  data.instr.Es.data(id,:) = inv(A) * data.instr.Es.data(id,:)';
endfor
#----Correction for Li
A = I +D_Li;
disp('** Straylight correction for Li,wait!')
for id = 1:length(data.instr.Li.data(:,1));
  data.instr.Li.data(id,:) = inv(A) * data.instr.Li.data(id,:)';
endfor
#----Correction for Lt
A = I +D_Lt;
disp('** Straylight correction for Lt,wait!')
for id = 1:length(data.instr.Lt.data(:,1));
  data.instr.Lt.data(id,:) = inv(A) * data.instr.Lt.data(id,:)';
endfor


endfunction


