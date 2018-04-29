# polygon-drawing
ecs175 graphics polygon drawing and coloring

This program takes coordinates from the data file "test_scene," and it will use DDA/Bresenham line algorithms to draw polygons from them. I tried to implement scan line rasterization, but it's not complete. Then the program lets you select polygons and do transformations with them like scaling, rotation, and translation. 

It'll ask you for: 
  grid specifications               (how big you want the window to be), 
  viewport (xmin,ymin,xmax,ymax)    (a square of area that's actually viewable to you), 
  scale factor, 
  rotation angle                    (how much to rotate by ever tick), 
  translation x and y.

Select a Polygon:
In "test_scene," the first line has the number of polygons, and the order the polygons are listed in this data file is the # ID of each polygon.
eg. If the top number is 5, then to target a polygon, you could press a number 0-4 on the keyboard.

Rotation: (Press A or S)
A for go degrees left
S for degrees right      

Scaling: (I did scaling by same factors on x and y)
W for increase by 120% default
Z for decrease by 80% default

Translation:
I,J,K,L -> Up, Left, Down, Right (It's like WASD)

Quit:
Q

Switch Bresenham and DDA:
P (will tell you which algorithm it's on once you press p)

Switch color or no color:
C    (see below)

The coloring is not complete. It's a bug that should be fixed. Rotation needs to be fixed a little, too because when I rotate, it moves a little left/right rather than just rotating
