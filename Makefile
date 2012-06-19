SOURCES = Main.cpp HodgkinHuxley.cpp
INCLUDES = HodgkinHuxley.hpp

Main: $(SOURCES) $(INCLUDES)
	g++ -o Main $(SOURCES) 