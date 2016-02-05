using Metrics, KMedians

function main()
	x1 = [-100 -100 ; -99 -99 ; -98 -98 ; -101 -101 ; -100 -101]
	x2 = -x1
	x = [x1 ; x2 ; (x1 - 5) ; (x2 + 5)]

	#sol = exact_kmedian(x', 4, 2, 1, 20, verbose = true)
	sol = kmedian_round(x', 4, 2, 1, 20, verbose = true)

	println("---------\nCENTERS\n-----------")
	println(sol.centers)
	println("---------\nASSIGNMENTS\n-----------")
	for x in sol.assignments ; println("\t  $(x)") ; end
	println("---------\nCLUSTERS\n-----------")
	for c in keys(sol.centers) ; println("$(c):\t$(sol.clusters[c])")  ; end

	

end


main()

