# PyRagdoll2.6_PyODE_PyOpenGL

1. sudo apt-get install freeglut3-dev # instala OpenGL necesario para compilar demos

make:
  g++ demo_basket.cpp -L/home/david/ode-0.16.2/ode/src/libode.la -L/home/david/ode-0.16.2/drawstuff/src/libdrawstuff.la -L/usr/lib -L/usr/local/lib -I/usr     /include -I/home/david/ode-0.16.2/include -lode -lGL -lGLU -lm
