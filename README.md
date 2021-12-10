Based on nialltl's (https://twitter.com/nialltl) 2D implementation of MLS-MPM: (https://github.com/nialltl/incremental_mpm).

Expanded to 3D and added liquid rendering techniques based on this Nvidia paper: https://developer.download.nvidia.com/presentations/2010/gdc/Direct3D_Effects.pdf

Current state:

https://www.youtube.com/watch?v=6VQqMkqwCZU&ab_channel=JackBlazes

Main current issue is with how Unity uses depth textures - I need to get one and apply a blur to it. At first I tried all sorts of methods to get Unity's built in depth texture representation and sample it - unfortunately this did not work. Similarly, when I tried manually creating a depth texture using a separate shader and camera, I ran into the same issue - it would not write to the new texture. This is because of the way the spheres are indirectly instanced. Would appreciate tips on any possible solutions.

6.839 Final Project
