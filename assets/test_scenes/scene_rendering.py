import bpy

X = bpy.data.objects["x"]
Y = bpy.data.objects["y"]
Z = bpy.data.objects["z"]

axis = [X, Y, Z]
scene = bpy.context.scene
cameras = [
    bpy.data.objects["Camera"],
    bpy.data.objects["Camera2"]
]

for ax in axis:
    ax.hide_render = True

for (i, ax) in enumerate(axis):
    ax.hide_render = False
    
    for (j, camera) in enumerate(cameras):
        scene.render.filepath = f"//renders/{j}/{i}.png"
        bpy.context.scene.camera = camera
        bpy.ops.render.render(write_still=True)
    
    ax.hide_render = True

for ax in axis:
    ax.hide_render = False