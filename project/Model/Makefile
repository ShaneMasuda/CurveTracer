CFLAGS= -Ofast -std=c++11
SRC= Metric.cpp LightPath.cpp ../CurveTracer/Polynomial.cpp ../CurveTracer/Primitives.cpp $(wildcard alglib/*.cpp)
HEADER= Metric.h Timer.h LightPath.h ../CurveTracer/Polynomial.h ../CurveTracer/Primitives.h Vector.h $(wildcard alglib/*.h)

all: sample test


sample: $(SRC) $(HEADER) sample.cpp
	g++ sample.cpp $(SRC) -o sample $(CFLAGS)

test: $(SRC) $(HEADER) test.cpp
	g++ test.cpp $(SRC) -o test $(CFLAGS)

messingaround: $(SRC) $(HEADER) messingaround.cpp
	g++ messingaround.cpp $(SRC) -o messingaround $(CFLAGS)

clean:
	rm *~ sample test
