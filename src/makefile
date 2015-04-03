optimization = -O3
objects = CmdLine.o PhenoBF.o CountData.o PhyloPop_params.o gzstream.o
libs = -L/opt/local/lib -lgsl -lgslcblas -lm -lz
inc = -I/opt/local/include


pheno: PhenoSelect.o $(objects)
	g++ $(libs) -o phenoselect PhenoSelect.o $(objects) $(optimization)

test: test.o $(objects)
	g++ $(libs) -o test test.o $(objects) $(optimization)

%.o: %.cpp
	g++ -c $< -o $@ $(inc) $(optimization)

clean:
	rm *.o
