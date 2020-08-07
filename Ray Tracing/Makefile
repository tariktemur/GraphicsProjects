CXX=g++
CXXFLAGS=-Wall -Wextra -pedantic -std=c++17 -pthread -O3

.PHONY: all
all: raytracer

.PHONY: clean
clean:
	rm *.o *.gch

Camera.o: Camera.cpp Camera.h Ray.h VectorUtility.hpp defs.h
	$(CXX) $(CXXFLAGS) -c Camera.cpp

Image.o: Image.cpp Image.h defs.h
	$(CXX) $(CXXFLAGS) -c Image.cpp

Light.o: Light.cpp Light.h defs.h
	$(CXX) $(CXXFLAGS) -c Light.cpp

main.o: main.cpp Scene.h Ray.h defs.h
	$(CXX) $(CXXFLAGS) -c main.cpp

Material.o: Material.cpp Material.h defs.h
	$(CXX) $(CXXFLAGS) -c Material.cpp

Ray.o: Ray.cpp Ray.h VectorUtility.hpp defs.h
	$(CXX) $(CXXFLAGS) -c Ray.cpp

Scene.o: Scene.cpp Scene.h Ray.h Camera.h Light.h Material.h Shape.h Image.h defs.h
	$(CXX) $(CXXFLAGS) -c Scene.cpp

Shape.o: Shape.cpp Shape.h Scene.h Ray.h VectorUtility.hpp defs.h
	$(CXX) $(CXXFLAGS) -c Shape.cpp

testscript.o: testscript.cpp Ray.h VectorUtility.hpp defs.h
	$(CXX) $(CXXFLAGS) -c testscript.cpp

testscript: testscript.o Ray.o VectorUtility.o
	$(CXX) -o testscript testscript.o Ray.o VectorUtility.o

tinyxml2.o: tinyxml2.cpp tinyxml2.h
	$(CXX) -std=c++17 -O3 -c tinyxml2.cpp

VectorUtility.o: VectorUtility.cpp VectorUtility.hpp defs.h
	$(CXX) $(CXXFLAGS) -c VectorUtility.cpp

raytracer: Camera.o Image.o Light.o main.o Material.o Ray.o Scene.o Shape.o tinyxml2.o VectorUtility.o
	$(CXX) -pthread -o raytracer Camera.o Image.o Light.o main.o Material.o Ray.o Scene.o Shape.o tinyxml2.o VectorUtility.o

.PHONY: render
render:
	./render.sh

.PHONY: test
test: testscript
	./testscript
