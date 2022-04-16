CPPFLAGS    = -I./inc -I./3rdParty/armadillo-10.1.2/include 
CPPFLAGS    += -I./3rdParty/knn-cpp/include 
CPPFLAGS    += -I/usr/include/eigen3
CPPFLAGS    += -I/usr/include/opencv4
CPPFLAGS    += -I./3rdParty/pugixml/src
LDFLAGS	    = -L/usr/lib/x86_64-linux-gnu/
LDFLAGS	    += -lopencv_core -lopencv_imgcodecs 
LDFLAGS     += -lboost_program_options 
LDFLAGS     += -L./3rdParty/pugixml/build -lpugixml
CXX         := g++ -O3 --std=c++11 


all: create wm ig aptdec

wm: 
	$(CXX) $(CPPFLAGS) src/web_merc.cpp -o bin/wm $^ $(LDFLAGS)

ig:
	$(CXX) $(CPPFLAGS) src/im_grid.cpp -o bin/ig $^ $(LDFLAGS)

aptdec: ; cd 3rdParty/aptdec && make && cp aptdec ../../bin

create: ; mkdir -p bin

clean: ; rm -rf bin
