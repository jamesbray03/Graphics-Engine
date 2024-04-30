# Graphics Engine

![high_res](https://github.com/jamesbray03/Graphics-Engine/assets/47334864/9e81e6f1-5d54-43c3-838a-3746ea5350bc)

A simple 3D graphics engine using [LML](https://github.com/jamesbray03/Lightweight-Matrix-Library) which renders obj files. 

## Setting up a scene

### Settings

Export settings:
- `render name`: The name of the file to render to (exclude folders and file extension).
- `screen width`: The width of the screen in pixels.
- `screen height`: The height of the screen in pixels.

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


### Example
```
--------------- SCENE SETUP ---------------

Export Settings
        render name: rendered_sphere
        screen width: 800
        screen height: 800

Camera Settings
    use default camera? (y/n): n
        perspective: 1
        near clipping plane: 0.1
        far clipping plane: 10
        field of view: 60

Object Settings 
    Object 1:   
        object: Sphere
        use default position? (y/n): n
            scale: 1 1 1
            rotation: 0 0 0
            translation: 0 0 3
        add another object? (y/n): n

--------------- FINISHED SETUP ---------------

Press any key to render...
```
