# Graphics Engine


A simple 3D graphics engine using [LML](https://github.com/jamesbray03/Lightweight-Matrix-Library) which renders obj files. 

## Setting up a scene

Export settings:
- `render name`: The name of the file to render to (exclude folders and file extension).
- `screen width`: The width of the screen in pixels. For ASCII, use about 100. For PGM, use whatever image resolution you like.
- `screen height`: The height of the screen in pixels. For ASCII, use about 80. For PGM, use whatever image resolution you like.

Camera settings:
- `perspective`: Set to 0 for orthographic projection, or 1 for perspective projection.
- `far clipping plane`: The distance from the camera to the far clipping plane (default: 10).
- `near clipping plane`: The distance from the camera to the near clipping plane (default: 0.1).
- `field of view`: The field of view in degrees (only used for perspective projection).

Object settings:
- `obj file`: The name of the obj file in the obj folder (exclude folders and file extension).
- `scale`: The scale of the object in world space (default: 1 1 1).
- `rotation`: The rotation of the object in world space (default: 0 0 0).
- `position`: The position of the object in world space (default: 0 0 8).
    Note: The z-axis is forward and the y-axis is up.


