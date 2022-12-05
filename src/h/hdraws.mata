capture mata: mata drop hdraws()
mata matrix function hdraws(N,D,P,B,A,I) {
	if (A==1) D 		= D/2 ;
	hdrawsvec 			= J(N*D,1,.)
	hdrawsvecburn		= ghalton(N*D+B,P,0)
	for (i=1; i<=N*D; i++) {
		hdrawsvec[i,1]	= hdrawsvecburn[i+B,1]
	}
	_hdraws				= colshape(hdrawsvec,D)
	if (A==1) hdraws	= _hdraws,(J(N,D,1) - _hdraws)
	else hdraws			= _hdraws
	names				= st_tempname(cols(hdraws))
	for (i=1; i<=cols(hdraws); i++) {
		st_addvar("double", names[i])
		st_store(.,names[i], hdraws[,i])
	}
	st_matrix(I,st_varindex(names))
}
mata: mata mlib create lrfrontier, dir("`c(sysdir_plus)'l/") replace
mata: mata mlib add lrfrontier hdraws(), dir("`c(sysdir_plus)'l/") complete