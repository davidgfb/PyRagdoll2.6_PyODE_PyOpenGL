# PyRagdoll2.6_PyODE_PyOpenGL

0. https://bitbucket.org/odedevs/ode/downloads/ode-0.16.2.tar.gz 

1. sudo apt-get install freeglut3-dev # instala OpenGL necesario para compilar demos # ./configure --prefix=/usr/local/

2. make

#############################################################################################################
g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_buggy.o -MD -MP -MF .deps/demo_buggy.Tpo -c -o demo_buggy.o demo_buggy.cpp
mv -f .deps/demo_buggy.Tpo .deps/demo_buggy.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link g++  -g -O2     -o demo_buggy demo_buggy.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread
libtool: link: g++ -g -O2 -o demo_buggy demo_buggy.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread




#############################################################################################################








make:
  g++ -g -O2 -o demo_trimesh demo_trimesh.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread
  
  g++  -g -O2     -o demo_trimesh demo_trimesh.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread

g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_trimesh.o -MD -MP -MF .deps/demo_trimesh.Tpo -c -o demo_trimesh.o demo_trimesh.cpp
mv -f .deps/demo_trimesh.Tpo .deps/demo_trimesh.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link

g++ -g -O2 -o demo_moving_convex demo_moving_convex.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread

g++  -g -O2     -o demo_moving_convex demo_moving_convex.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread

g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_moving_convex.o -MD -MP -MF .deps/demo_moving_convex.Tpo -c -o demo_moving_convex.o demo_moving_convex.cpp
mv -f .deps/demo_moving_convex.Tpo .deps/demo_moving_convex.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link 

g++ -g -O2 -o demo_moving_trimesh demo_moving_trimesh.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread

g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_moving_trimesh.o -MD -MP -MF .deps/demo_moving_trimesh.Tpo -c -o demo_moving_trimesh.o demo_moving_trimesh.cpp
mv -f .deps/demo_moving_trimesh.Tpo .deps/demo_moving_trimesh.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link g++  -g -O2     -o demo_moving_trimesh demo_moving_trimesh.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread

g++ -g -O2 -o demo_cyl demo_cyl.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread

g++  -g -O2     -o demo_cyl demo_cyl.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread

g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_cyl.o -MD -MP -MF .deps/demo_cyl.Tpo -c -o demo_cyl.o demo_cyl.cpp
mv -f .deps/demo_cyl.Tpo .deps/demo_cyl.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link

# lo de arriba esta copiado de abajo a arriba. Dale la vuelta

####### demo basket #############

g++ -DHAVE_CONFIG_H -I. -I../../ode/src  -I../../include -I../../include -DDRAWSTUFF_TEXTURE_PATH="\"/home/david/ode-0.16.2/drawstuff/textures\"" -DdTRIMESH_ENABLED   -g -O2 -MT demo_basket.o -MD -MP -MF .deps/demo_basket.Tpo -c -o demo_basket.o demo_basket.cpp
mv -f .deps/demo_basket.Tpo .deps/demo_basket.Po
/bin/bash ../../libtool  --tag=CXX   --mode=link

g++  -g -O2     -o demo_basket demo_basket.o ../../drawstuff/src/libdrawstuff.la ../../ode/src/libode.la -lGLU -lGL  -lrt -lm  -lpthread

g++ -g -O2 -o demo_basket demo_basket.o  ../../drawstuff/src/.libs/libdrawstuff.a -lX11 ../../ode/src/.libs/libode.a -lGLU -lGL -lrt -lm -lpthread

####### fin demo basket #########

  #g++ demo_basket.cpp -L/home/david/ode-0.16.2/ode/src/libode.la -L/home/david/ode-0.16.2/drawstuff/src/libdrawstuff.la -L/usr/lib -L/usr/local/lib -I/usr   #  /include -I/home/david/ode-0.16.2/include -lode -lGL -lGLU -lm
