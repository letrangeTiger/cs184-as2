CC = g++
ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
  CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
  LDFLAGS = -framework GLUT -framework OpenGL \
	        -L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
	        -lGL -lGLU -lm -lstdc++
  else
  CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
  LDFLAGS = -lglut -lGLU
endif

RM = /bin/rm -f
all: main
main: Scene.o
	$(CC) $(CFLAGS) -o as2 Scene.o $(LDFLAGS) 
Scene.o: Scene.cpp
	$(CC) $(CFLAGS) -c Scene.cpp -o Scene.o
clean:
	$(RM) *.o as2
