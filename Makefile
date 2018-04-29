CC = /usr/bin/g++

#OPENGL LIBS FOR LINUX
GLLIB := -L/usr/lib/ -lGL -lGLU -lGLEW -lglut
#OPENGL LIBS FOR MAC
#GLLIB := -framework OpenGL -framework GLUT

#COMPILER FLAGS
CCFLAGS :=

#include directories
#should include gl.h glut.h etc...
INCDIR := -I/usr/include
LDLIBS := $(GLLIB)

TARGET = drawing
OBJS = drawing.o


all: $(TARGET)


$(TARGET): $(OBJS)
	$(CC) $^ $(CCFLAGS) $(LDLIBS) -o $@

%.o : %.cpp
	$(CC) $(CCFLAGS) -o $@ -c $(INCDIR) $<

clean:
	rm -rf $(OBJS) $(TARGET)

