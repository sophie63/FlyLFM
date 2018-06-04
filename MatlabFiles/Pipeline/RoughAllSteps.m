function RoughAllSteps(file)
%file= '/media/test/100628regc.nii';
Dkf= modified_dFF_psf_KF(file);
Dica = A_ICAalaMelodicfunc(file,Dkf,100);
Dcolor = ColoredSumaryMapfunc(file,Dica);
