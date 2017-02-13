
diamond:
	diamond -s /data/PhD/fluidity_trunk/schemas/fluidity_options.rng lock_exchange_2d_cg.flml

input: clean
	gmsh -2 ./channel_2d.geo
	#cp src/channel_2d.msh .
	dolfin-convert channel_2d.msh channel_2d.xml

run: input
	fluidity lock_exchange_2d_cg.flml

clean: 
	rm -f  ./src/*.msh *.msh channel_2d*.xml *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump*
