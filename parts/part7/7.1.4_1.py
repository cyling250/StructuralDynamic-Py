# 导入依赖库
from ast import mod
from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

# 设定材料参数
beam_width = 0.3  # 梁宽度为0.3m
beam_depth = 0.6  # 梁厚度为0.6m
span_length = 12  # 梁跨度为12m
E = 209e9  # 弹性模量，单位：Pa (Pascal)
mu = 0.3  # 泊松比
load_magnitude = 12000.0  # 荷载大小，单位：N (牛顿)
element_size = 0.1  # 网格单元尺寸，单位：m
# 创建新模型和新视口
my_model = mdb.Model(name='SimpleSupportedBeam')
my_viewport = session.Viewport(name='SimpleSupportedBeam_View')
my_viewport.makeCurrent()
# 创建基础特征草图
my_sketch = my_model.ConstrainedSketch(name='BeamSketch', sheetSize=10.0)
# 绘制梁截面
my_sketch.Line(point1=(0.0, 0.0), point2=(beam_width, 0.0))  # 底边
my_sketch.Line(point1=(beam_width, 0.0), point2=(beam_width, beam_depth))  # 右边
my_sketch.Line(point1=(beam_width, beam_depth), point2=(0.0, beam_depth))  # 顶边
my_sketch.Line(point1=(0.0, beam_depth), point2=(0.0, 0.0))  # 左边
# 创建三维变形体部件
my_part = my_model.Part(name='SimpleSupportedBeam_Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
my_part.BaseSolidExtrude(sketch=my_sketch, depth=span_length)  # 对草图拉伸创建部件
# 创建材料
my_material = my_model.Material(name='Steel')  # 创建材料为钢材
elastic_properties = (E, mu)
my_material.Elastic(table=(elastic_properties,))
# 创建实体截面
my_section = my_model.HomogeneousSolidSection(name='BeamSection', material='Steel', thickness=None)
# 创建简支梁截面
region = (my_part.cells,)
my_part.SectionAssignment(region=region, sectionName='BeamSection', offset=0.0,
                          offsetType=MIDDLE_SURFACE, offsetField='',
                          thicknessAssignment=FROM_SECTION)
# 创建装配
my_assembly = my_model.rootAssembly
my_instance = my_assembly.Instance(name='SimpleSupportedBeam_Instance', part=my_part, dependent=ON)
# 创建分析步
my_model.StaticStep(name='beam_load', previous='Initial')
# 找到荷载施加的面
region_load_center = (0.1, 0.6, 6)
region_load = my_instance.faces.findAt((region_load_center,))
region_load = my_assembly.Surface(side1Faces=region_load, name="load_surface")
# 找到约束施加的面
region_fix_center = (0.1, 0.6, 0)
region_fix = my_instance.faces.findAt((region_fix_center,))
region_fix = my_assembly.Set(faces=region_fix, name="fix_Set")
# 创建均布荷载
my_model.Pressure(name='Load',
                  createStepName='beam_load', region=region_load, distributionType=UNIFORM,
                  field='', magnitude=load_magnitude, amplitude=UNSET)
# 创建固支约束
my_model.EncastreBC(name='FixedSupport', createStepName='beam_load', region=region_fix)
# 创建网格
my_part.seedPart(size=element_size)
my_part.generateMesh()
# 提交作业
my_job = mdb.Job(name='SimpleSupportedBeam_Job', model='SimpleSupportedBeam',
                 description='Simple Supported Beam Analysis')
my_job.submit(consistencyChecking=OFF)
mdb.jobs[my_job.name].waitForCompletion()
# 汇总数据，打开OBD文件
session.mdbData.summary()
session.openOdb(name='F:/temp/SimpleSupportedBeam_Job.odb')
