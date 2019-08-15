
# Camera and film

A "camera" is a system that shoots _vision rays_ at a scene. These rays are traced until
they hit the surface of an object (or miss the object and escape into infinity). The
contributions of all light sources to surface brightness at that location are then added up
and added to the "film" associated with the camera.

The film keeps track of surface brightness "seen" by the camera in some way. The film is,
in a sense, what is "returned" by the whole rendering process. The contents of the film
can then be saved to disk in some way.

* `ImageFilm` keeps a pixel grid and creates a JPEG image of the scene as viewed
  through the camera.
* `PhotometricFilm` simply adds up all the energy received by the camera, being effectively
  a single-pixel detector.
* `DelayFilm` notes the distances travelled by each ray, and keeps a histogram of
  brightness as a function of distance into the scene. It is used to simulate time-resolved
  observations with laser ranging


## Orthographic camera

This is the simplest camera to use. It views the scene in orthographic projection. It is
simply a rectangular "window" in the scene, shooting out parallel vision rays.

To create a camera and a film to make a JPEG image, we define the pixel resolution of the
film as well as the size of the window. Then we need to rotate and shift the camera to
point it where we want.

Note that the resolution is a property of the film and the window is a property of the
camera. These should match so that the ratio of the X and Y resolutions should match the
ratio of the X and Y widths of the window.

```
resX = 1024
resY = 768

F = ImageFilm(resX, resY, 1.0, TriangleFilter(1, 1))

widthX = 20.0
widthY = 15.0

camera_position = rotation(Y_AXIS, angle) * translation(0.0, 0.0, -25.0)

image_cam = OrthographicCamera(camera_position, widthX, widthY, 0.0, 0.0, F)
```

By default, the camera is located at the origin and is looking in the +Z direction. To view
an object at the origin, the camera is first shifted 25 metres in the -Z direction and then
rotated around the Y axis.

TODO: A utility function to generate the transformation to place a camera at a given
point, looking at another given point.
