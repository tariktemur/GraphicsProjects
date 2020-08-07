CXX=g++
CXXFLAGS=-Wall -Wextra -pedantic -std=c++17 -O3

.PHONY: all
all: rasterizer

.PHONY: clean
clean:
	rm *.o *.gch *.ppm vgcore.*

Camera.o: Camera.cpp Camera.h Vec3.h
	$(CXX) $(CXXFLAGS) -c $<

Color.o: Color.cpp Color.h
	$(CXX) $(CXXFLAGS) -c $<

Helpers.o: Helpers.cpp Helpers.h Matrix4.h Vec3.h Vec4.h
	$(CXX) $(CXXFLAGS) -c $<

Main.o: Main.cpp Scene.h Matrix4.h Helpers.h
	$(CXX) $(CXXFLAGS) -c $<

Matrix4.o: Matrix4.cpp Matrix4.h
	$(CXX) $(CXXFLAGS) -c $<

Model.o: Model.cpp Model.h Triangle.h
	$(CXX) $(CXXFLAGS) -c $<

Rotation.o: Rotation.cpp Rotation.h
	$(CXX) $(CXXFLAGS) -c $<

Scaling.o: Scaling.cpp Scaling.h
	$(CXX) $(CXXFLAGS) -c $<

Scene.o: Scene.cpp Scene.h Camera.h Color.h Matrix4.h Model.h Rotation.h Scaling.h Translation.h Triangle.h Vec3.h tinyxml2.h Helpers.h
	$(CXX) $(CXXFLAGS) -c $<

tinyxml2.o: tinyxml2.cpp tinyxml2.h
	$(CXX) -std=c++17 -O3 -c $<

Translation.o: Translation.cpp Translation.h
	$(CXX) $(CXXFLAGS) -c $<

Triangle.o: Triangle.cpp Triangle.h
	$(CXX) $(CXXFLAGS) -c $<

Vec3.o: Vec3.cpp Vec3.h Vec4.h
	$(CXX) $(CXXFLAGS) -c $<

Vec4.o: Vec4.cpp Vec4.h Vec3.h
	$(CXX) $(CXXFLAGS) -c $<

rasterizer: Camera.o Color.o Helpers.o Main.o Matrix4.o Model.o Rotation.o Scaling.o Scene.o tinyxml2.o Translation.o Triangle.o Vec3.o Vec4.o
	$(CXX) -O3 -flto -o $@ $^
