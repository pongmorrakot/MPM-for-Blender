import bpy
#import mpm
from bpy import context
from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       EnumProperty,
                       )

class MPM_BakeButton(bpy.types.Operator):
    bl_idname = "wm.mpm_bake"
    bl_label = "Bake"

    def execute(self, context):
        print("Bake")
        # call the main mpm function which will then call the parser that will parse all the data
        mpm.parse(context.selected_objects,
                bpy.data.objects['domain'],
				[context.scene.my_addon.physics_x,context.scene.my_addon.physics_y, context.scene.my_addon.physics_z],
                context.scene.my_addon.gridRes,
                context.scene.my_addon.young,
                context.scene.my_addon.compression,
                context.scene.my_addon.stretch,
                context.scene.my_addon.hardening,
				context.scene.my_addon.simLength,
				context.scene.my_addon.substep)
        return {'FINISHED'}

class MPM_Property(bpy.types.PropertyGroup):
    bl_idname = "PARTICLE_OT_mpm_property"
    bl_label = "MPM Property"
    bl_options = {'REGISTER', 'UNDO'}
    
    enable: bpy.props.BoolProperty(
        name="Enable",
        description="",
        default = False)
    
    physics_x: bpy.props.FloatProperty(
        name="X",
        description="",
        default = 0)
    physics_y: bpy.props.FloatProperty(
        name="Y",
        description="",
        default = -1)
    physics_z: bpy.props.FloatProperty(
        name="Z",
        description="",
        default = 0)
        
    young: bpy.props.IntProperty(
        name="Initial Young Modulus",
        description="",
        default = 5000,
        min = 0)
    compression: bpy.props.FloatProperty(
        name="Critical Compression",
        description="",
        default = 0.025,
        min = 0)
    stretch: bpy.props.FloatProperty(
        name="Critical Stretch",
        description="",
        default = 0.0075,
        min = 0)
    hardening: bpy.props.FloatProperty(
        name="Hardening Coefficient",
        description="",
        default = 10,
        min = 0)
        
    gridRes: bpy.props.IntProperty(
        name="Grid Resolution",
        description="",
        default = 10,
        min = 0)
    simLength: bpy.props.IntProperty(
        name="Simulation Length(second)",
        description="",
        default = 12,
        min = 0)
    substep: bpy.props.IntProperty(
        name="Substep",
        description="",
        default = 1,
        min = 1)
    


class MPM_Panel(bpy.types.Panel):
    bl_idname = "PARTICLE_PT_mpm_panel"
    bl_label = "Material Point Method"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "particle"
    

    def draw(self, context):
        # You can set the property values that should be used when the user
        # presses the button in the UI.
        layout = self.layout
        object = context.object
        
        row = layout.row()
        row.prop(context.scene.my_addon, "enable")
        
        if context.scene.my_addon.enable == True:
            layout.label(text="Material Point Method")
            layout.label(text="Physics")
            row = layout.row()
            row.prop(context.scene.my_addon, "physics_x")
            row.prop(context.scene.my_addon, "physics_y")
            row.prop(context.scene.my_addon, "physics_z")
            layout.label(text="Properties")
            row = layout.row()
            row.prop(context.scene.my_addon, "young")
            row = layout.row()
            row.prop(context.scene.my_addon, "compression")
            row = layout.row()
            row.prop(context.scene.my_addon, "stretch")
            row = layout.row()
            row.prop(context.scene.my_addon, "hardening")
            
            layout.label(text="Simulation")
            row = layout.row()
            row.prop(context.scene.my_addon, "gridRes")
            row = layout.row()
            row.prop(context.scene.my_addon, "simLength")
            row = layout.row()
            row.prop(context.scene.my_addon, "substep")
            layout.label(text="")
            row = layout.row()
            row.operator(MPM_BakeButton.bl_idname)


def register():
    bpy.utils.register_class(MPM_BakeButton)
    bpy.utils.register_class(MPM_Property)
    bpy.utils.register_class(MPM_Panel)
    bpy.types.Scene.my_addon = bpy.props.PointerProperty(type=MPM_Property)
        
def unregister():
    bpy.utils.unregister_class(MPM_BakeButton)
    bpy.utils.unregister_class(MPM_Property)
    bpy.utils.unregister_class(MPM_Panel)
    
if __name__ == "__main__" :  
    register()	