

function observabilidad(A,C)

	Obs = obsv(A,C);
	rango_obs = rank(Obs);
	cant_estados = rank(A);
	estados_no_observables = cant_estados - rango_obs

end


