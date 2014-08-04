%% OOFEM generated export file 
% Output for time 2.500000
function [mesh area data specials ReactionForces IntegrationPointFields]=tr2shell7XFEM_CZ_out_m0_6

	mesh=[];
	data=[];
	area.area=[];
	area.volume=[];
	specials=[];

 %% Export of reaction forces 

	ReactionForces.DofManNumbers = [1 2 3 4 5 6 ];
	ReactionForces.ReactionForces = cell(6,1); 
	ReactionForces.DofIDs = cell(6,1); 
	ReactionForces.ReactionForces{1} = [6.711511e-012 7.437079e-012 -1.692230e-014 1.662544e-015 1.645220e-015 -1.532163e-012 -1.394796e-017 3.355755e-015 3.718540e-015 -8.681382e-018 8.481063e-017 9.383659e-017 -7.664311e-016 ];
	ReactionForces.DofIDs{1} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	ReactionForces.ReactionForces{2} = [1.875748e-012 3.764932e-050 -1.111795e-014 -1.461538e-015 -1.476540e-016 2.212494e-013 -8.583992e-018 9.378738e-016 3.170113e-053 -7.214362e-018 2.187544e-017 2.646451e-019 1.104067e-016 ];
	ReactionForces.DofIDs{2} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	ReactionForces.ReactionForces{3} = [5.684372e-028 2.078531e-012 -3.027681e-015 4.402021e-016 -1.476540e-016 2.211621e-013 -2.962222e-018 1.387786e-031 1.039266e-015 -5.922499e-018 -6.802236e-021 2.470890e-017 1.105092e-016 ];
	ReactionForces.DofIDs{3} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	ReactionForces.ReactionForces{4} = [-8.587258e-012 1.201486e-012 2.212459e-014 -6.083033e-016 5.906160e-016 -8.854395e-013 1.517687e-017 -4.293629e-015 6.007431e-016 1.293085e-017 -1.078152e-016 1.482995e-017 -4.423239e-016 ];
	ReactionForces.DofIDs{4} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	ReactionForces.ReactionForces{5} = [-1.084268e-012 -1.201486e-012 3.595903e-015 -1.156447e-015 -1.087153e-015 -9.967430e-014 4.830663e-018 -5.421340e-016 -6.007431e-016 -2.096756e-018 -1.444195e-017 -1.573259e-017 -4.972554e-017 ];
	ReactionForces.DofIDs{5} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	ReactionForces.ReactionForces{6} = [1.084268e-012 -9.515611e-012 5.347443e-015 2.865839e-015 5.906160e-016 -8.851503e-013 5.486635e-018 5.421340e-016 -4.757805e-015 1.098415e-017 1.378053e-017 -1.193863e-016 -4.424431e-016 ];
	ReactionForces.DofIDs{6} = [1 2 3 15 16 17 18 23 24 25 26 27 28 ];
	IntegrationPointFields=[];

end
