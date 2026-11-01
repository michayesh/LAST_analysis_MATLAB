VisitPath = '/mnt/euclid/last/data/LAST.01.01.01/2026/07/22/proc/200041v0/';
% VisitPath = '/bigdata/projects/last/data/LAST.01.01.01/2026/07/22/proc/200041v0/';
% /bigdata/projects/last/data/LAST.01.01.01/2026/06

[~, Coadd] = pipeline.last.load.loadVisit(VisitPath, 'TempName_IndivIm',[], 'TempName_MS',[], 'Ncrop',0);
[Coadd, ~] = imProc.psf.populatePSF(Coadd, 'RePopulatePSF',true,'Method','new', 'WingsMethod','empirical');