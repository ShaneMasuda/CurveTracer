CFLAGS= -g -Ofast -std=c++11 `libpng-config --ldflags`
SRC= Polynomial.cpp Primitives.cpp ../Render/LightPath.cpp ../Render/Metric.cpp ../Simplex/SimplexNoise.cpp $(wildcard ../Model/alglib/*.cpp)
HEADER= Polynomial.h Primitives.h ../Render/LightPath.h ../Render/Metric.h ../Render/Vector.h Camera.h ../Simplex/SimplexNoise.h ../Render/Timer.h $(wildcard ../Model/alglib/*.h) $(wildcard ../png++/*.hpp)

all: CurveTracer

CurveTracer: $(SRC) $(HEADER) CurveTracer.cpp
	g++ CurveTracer.cpp $(SRC) -o CurveTracer $(CFLAGS)

clean:
	rm *~ CurveTracer
