import bpy
import os
import tempfile

# ----------------------------
# CONFIG
# ----------------------------
CAMERAS = ("Camera1", "Camera2")
BACKGROUNDS = ("Background1", "Background2")
RENDERS = ("Render1", "Render2")
RESULTS = ("Result1", "Result2")

RENDER_PERCENT = 25  # render at lower percent for speed

# ----------------------------
# INTERNAL STATE
# ----------------------------
_state = {
    "is_processing": False,
}


# ----------------------------
# HELPERS
# ----------------------------
def get_image(name):
    img = bpy.data.images.get(name)
    if img is None:
        raise RuntimeError(f"Image '{name}' not found.")
    return img


def ensure_image(name, width, height, alpha=True):
    img = bpy.data.images.get(name)
    if img is None:
        img = bpy.data.images.new(name=name, width=width, height=height, alpha=alpha, float_buffer=True)
    else:
        if img.size[0] != width or img.size[1] != height:
            img.scale(width, height)
    return img


def alpha_over(bg_pixels, fg_pixels):
    # Both are flat RGBA arrays in straight alpha
    n = len(bg_pixels)
    out = [0.0] * n
    for i in range(0, n, 4):
        br, bg, bb, ba = bg_pixels[i:i+4]
        fr, fg, fb, fa = fg_pixels[i:i+4]

        oa = fa + ba * (1.0 - fa)

        if oa > 1e-8:
            or_ = (fr * fa + br * ba * (1.0 - fa)) / oa
            og = (fg * fa + bg * ba * (1.0 - fa)) / oa
            ob = (fb * fa + bb * ba * (1.0 - fa)) / oa
        else:
            or_ = og = ob = 0.0

        out[i] = or_
        out[i + 1] = og
        out[i + 2] = ob
        out[i + 3] = oa

    return out


def scale_to_canvas_match_width(src_pixels, src_w, src_h, dst_w, dst_h):
    # Scale source to match destination width, keep aspect ratio for height.
    # Then place it on a transparent dst_w x dst_h canvas (centered vertically).
    if src_w <= 0 or src_h <= 0:
        raise RuntimeError(f"Invalid source render size: {src_w}x{src_h}")

    scaled_h = max(1, int(round(dst_w * (src_h / src_w))))

    tmp_img = bpy.data.images.new(
        name="__overlay_scale_tmp__",
        width=src_w,
        height=src_h,
        alpha=True,
        float_buffer=True,
    )
    try:
        tmp_img.pixels.foreach_set(src_pixels)
        tmp_img.update()
        tmp_img.scale(dst_w, scaled_h)
        scaled_pixels = tmp_img.pixels[:]
    finally:
        bpy.data.images.remove(tmp_img)

    canvas = [0.0] * (dst_w * dst_h * 4)

    if scaled_h <= dst_h:
        dst_y0 = (dst_h - scaled_h) // 2
        src_y0 = 0
        rows = scaled_h
    else:
        dst_y0 = 0
        src_y0 = (scaled_h - dst_h) // 2
        rows = dst_h

    row_stride = dst_w * 4
    for r in range(rows):
        src_row_start = (src_y0 + r) * row_stride
        src_row_end = src_row_start + row_stride
        dst_row_start = (dst_y0 + r) * row_stride
        dst_row_end = dst_row_start + row_stride
        canvas[dst_row_start:dst_row_end] = scaled_pixels[src_row_start:src_row_end]

    return canvas


def render_camera_to_image(scene, camera_obj, target_render_img, width, height):
    # Save old render settings
    old_camera = scene.camera
    old_x = scene.render.resolution_x
    old_y = scene.render.resolution_y
    old_pct = scene.render.resolution_percentage
    old_filepath = scene.render.filepath
    old_format = scene.render.image_settings.file_format
    old_color_mode = scene.render.image_settings.color_mode

    try:
        scene.camera = camera_obj
        scene.render.resolution_x = width
        scene.render.resolution_y = height
        scene.render.resolution_percentage = max(1, min(100, int(getattr(scene, "overlay_render_percent", RENDER_PERCENT))))

        bpy.ops.render.render(write_still=False, use_viewport=False)

        rr = bpy.data.images.get("Render Result")
        rr_pixels = rr.pixels[:] if rr is not None else []
        rr_w, rr_h = (rr.size[0], rr.size[1]) if rr is not None else (0, 0)

        # Some Blender setups expose an empty Render Result buffer.
        # Fallback to rendering to a temporary PNG and loading it.
        if len(rr_pixels) != rr_w * rr_h * 4 or rr_w == 0 or rr_h == 0:
            fd, tmp_path = tempfile.mkstemp(suffix=".png")
            os.close(fd)
            loaded = None
            try:
                scene.render.filepath = tmp_path
                scene.render.image_settings.file_format = "PNG"
                scene.render.image_settings.color_mode = "RGBA"
                bpy.ops.render.render(write_still=True, use_viewport=False)

                loaded = bpy.data.images.load(tmp_path, check_existing=False)
                rr_pixels = loaded.pixels[:]
                rr_w, rr_h = loaded.size[0], loaded.size[1]

                if len(rr_pixels) != rr_w * rr_h * 4 or rr_w == 0 or rr_h == 0:
                    raise RuntimeError(
                        f"Rendered pixel buffer has unexpected size {len(rr_pixels)} "
                        f"for dimensions {rr_w}x{rr_h}."
                    )
            finally:
                if loaded is not None:
                    bpy.data.images.remove(loaded)
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

        # Ensure target is right size
        if target_render_img.size[0] != width or target_render_img.size[1] != height:
            target_render_img.scale(width, height)

        mapped_pixels = scale_to_canvas_match_width(rr_pixels, rr_w, rr_h, width, height)
        target_render_img.pixels.foreach_set(mapped_pixels)
        target_render_img.update()

    finally:
        scene.camera = old_camera
        scene.render.resolution_x = old_x
        scene.render.resolution_y = old_y
        scene.render.resolution_percentage = old_pct
        scene.render.filepath = old_filepath
        scene.render.image_settings.file_format = old_format
        scene.render.image_settings.color_mode = old_color_mode


def update_all_overlays():
    scene = bpy.context.scene

    for cam_name, bg_name, rnd_name, res_name in zip(CAMERAS, BACKGROUNDS, RENDERS, RESULTS):
        cam = bpy.data.objects.get(cam_name)
        if cam is None or cam.type != 'CAMERA':
            print(f"[Overlay] Skipping: camera '{cam_name}' not found or not a camera.")
            continue

        bg_img = get_image(bg_name)
        w, h = bg_img.size[0], bg_img.size[1]

        render_img = ensure_image(rnd_name, w, h, alpha=True)
        result_img = ensure_image(res_name, w, h, alpha=True)

        render_camera_to_image(scene, cam, render_img, w, h)

        bg_pixels = list(bg_img.pixels[:])
        fg_pixels = list(render_img.pixels[:])

        out_pixels = alpha_over(bg_pixels, fg_pixels)
        result_img.pixels.foreach_set(out_pixels)
        result_img.update()

    # Redraw image editors so updates are visible
    wm = bpy.context.window_manager
    for win in wm.windows:
        for area in win.screen.areas:
            if area.type == 'IMAGE_EDITOR':
                area.tag_redraw()

    print("[Overlay] Updated Result1/Result2.")


# ----------------------------
# MANUAL UPDATE UI
# ----------------------------
class OVERLAY_OT_update_results(bpy.types.Operator):
    bl_idname = "overlay.update_results"
    bl_label = "Update Overlay Results"
    bl_description = "Render Camera1/Camera2 and composite over Background1/Background2"

    def execute(self, context):
        if _state["is_processing"]:
            self.report({'INFO'}, "Overlay update already running")
            return {'CANCELLED'}

        _state["is_processing"] = True
        try:
            update_all_overlays()
        except Exception as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        finally:
            _state["is_processing"] = False

        self.report({'INFO'}, "Overlay results updated")
        return {'FINISHED'}


class OVERLAY_PT_tools(bpy.types.Panel):
    bl_label = "Line Overlay"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Overlay'

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        layout.prop(scene, "overlay_render_percent")
        layout.operator("overlay.update_results", icon='RENDER_STILL')


_CLASSES = (
    OVERLAY_OT_update_results,
    OVERLAY_PT_tools,
)


def cleanup_legacy_move_handlers():
    # If older versions of this script were run, remove their depsgraph handler(s)
    # so moving objects no longer triggers automatic renders.
    handlers = bpy.app.handlers.depsgraph_update_post
    to_remove = []
    for h in handlers:
        if getattr(h, "__name__", "") == "depsgraph_move_handler":
            to_remove.append(h)

    for h in to_remove:
        handlers.remove(h)

    if to_remove:
        print(f"[Overlay] Removed {len(to_remove)} legacy move handler(s).")


def unregister_ui():
    for cls in reversed(_CLASSES):
        try:
            bpy.utils.unregister_class(cls)
        except Exception:
            pass

    if hasattr(bpy.types.Scene, "overlay_render_percent"):
        del bpy.types.Scene.overlay_render_percent


def register_ui():
    cleanup_legacy_move_handlers()
    unregister_ui()

    bpy.types.Scene.overlay_render_percent = bpy.props.IntProperty(
        name="Render %",
        description="Temporary render scale before upscaling to texture width",
        default=RENDER_PERCENT,
        min=1,
        max=100,
    )

    for cls in _CLASSES:
        bpy.utils.register_class(cls)

    print("[Overlay] UI registered. Use N-panel > Overlay > Line Overlay > Update Overlay Results.")


# Backward compatibility with earlier script commands.
def unregister_handler():
    unregister_ui()


def register_handler():
    register_ui()


# ----------------------------
# RUN
# ----------------------------
register_ui()
