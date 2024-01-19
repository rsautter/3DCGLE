# Stochastic Complex Ginzburg-Landau

Implementation of the Complex Ginzburg-Landau Equation (CGLE)  using pseudospectral method and Runge-Kuta-Fehlberg 4-5 method.
The additive CGLE is:
    $$\partial_t A = (1+ib) \nabla^2 A + A  - (1+ic) |A|^2A$$
where $A$ is a complex number.

The following video shows the  traditional Complex Ginzburg-Landau:

https://user-images.githubusercontent.com/14216783/224129416-69f13958-8244-4741-9cf3-d6b9e51144e2.mp4

The following video shows the 3D amplitudes of a CGLE solution:

https://user-images.githubusercontent.com/14216783/224053529-af028b73-6f58-47f4-a0a4-0108392fa167.mp4

The following video shows the 3D gradient of the corresponding CGLE solution:

https://user-images.githubusercontent.com/14216783/224053529-af028b73-6f58-47f4-a0a4-0108392fa167.mp4

## GPA

We also present a Gradient Pattern Analysis ([GPA](https://en.wikipedia.org/wiki/Gradient_pattern_analysis)) of the system. The implementation is public available [here](https://github.com/rsautter/GPA). The 3D version of GPA is in the GPA3D folder.


## Paper link
Under submission
