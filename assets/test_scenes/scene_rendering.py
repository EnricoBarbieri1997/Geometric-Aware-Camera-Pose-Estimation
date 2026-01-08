import bpy

X = bpy.data.objects["x"]
Y = bpy.data.objects["y"]
Z = bpy.data.objects["z"]

axis = [X, Y, Z]
scene = bpy.context.scene

for ax in axis:
    ax.hide_render = True

for (i, ax) in enumerate(axis):
    ax.hide_render = False
    
    scene.render.filepath = f"//renders/{i}.png"
    bpy.ops.render.render(write_still=True)
    
    ax.hide_render = True

for ax in axis:
    ax.hide_render = False