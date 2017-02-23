msh: clean
	gmsh -2 unitsquaremesh.geo

xml: msh
	dolfin-convert unitsquaremesh.msh unitsquaremesh.xml

clean:
	rm -f  *.xml
