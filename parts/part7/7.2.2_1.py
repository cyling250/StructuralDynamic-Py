from ansys.mapdl import core as pymapdl
from ansys.mapdl.core import launch_mapdl
import numpy as np
print("No error means installing successfully.")
new_path = 'D:/ANASYS Inc/v231/ansys/bin/winx64/ANSYS231.exe'  # 请将这里改成你的ANSYS安装路径
pymapdl.change_default_ansys_path(new_path)

mapdl = launch_mapdl()
# create a rectangle with a few holes
mapdl.prep7()
rect_anum = mapdl.blc4(width=1, height=0.2)
# create several circles in the middle in the rectangle
for x in np.linspace(0.1, 0.9, 8):
    mapdl.cyl4(x, 0.1, 0.025)
# Generate a line plot
mapdl.lplot(color_lines=True, cpos='xy')

# -----------------------------------------------

plate_holes = mapdl.asba(rect_anum, 'all')
# extrude this area
mapdl.vext(plate_holes, dz=0.1)
mapdl.vplot()

# -----------------------------------------------

mapdl.et(1, 'SOLID186')
mapdl.vsweep('ALL')
mapdl.esize(0.1)
mapdl.eplot()
input("Press enter to exit.")
mapdl.exit()
