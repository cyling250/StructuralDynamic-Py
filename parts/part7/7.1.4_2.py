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

"""
本程序的相关库代码可以在libs文件夹中的abaqus.py中，为了在ABAQUS中调试方便，本文件也给出的了相关的函数。
需要注意，本文件的函数与abaqus.py中的代码是一致的
"""


def create_model(name="SimpleSupportBeam"):
    """
    创建模型
    Parameters
    ----------
    name 模型名称

    Returns 创建的模型
    -------

    """
    my_model = mdb.Model(name=name)
    # my_viewport = session.Viewport(name='SimpleSupportedBeam_View')
    # my_viewport.makeCurrent()
    return my_model


def beam_part(model, beam_width, beam_depth, length, name="beam_part", E=30e9, mu=0.2, material=None, section=None):
    """
    创建梁部件，并且该部件被分配材料与截面
    Parameters
    ----------
    model 模型名称
    beam_width 梁宽
    beam_depth 梁厚
    length 梁跨度
    name 梁部件名称
    E 弹性模量
    mu 泊松比
    material 材料名称
    section 截面名称

    Returns 创建的梁部件
    -------

    """
    # 创建草图
    sketch = model.ConstrainedSketch(name='beam_sketch', sheetSize=2 * max(beam_width, beam_depth))
    sketch.setPrimaryObject(option=STANDALONE)
    # 绘制草图
    sketch.rectangle(point1=(-beam_width / 2, -beam_depth / 2), point2=(beam_width / 2, beam_depth / 2))
    # 创建部件
    part = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=sketch, depth=length)
    # 如果没有指定该梁的材料或者材料没有被定义，则自定义材料；当材料被定义时，弹性模量E与泊松比mu没用
    if not material:
        material = model.Material(name='concrete')
        material.Elastic(table=((E, mu),))
    # 如果没有指定该梁的截面或者截面为被定义，则自定义截面
    if not section:
        section = model.HomogeneousSolidSection(name='beam_section', material='concrete', thickness=None)
    # 分配截面到实体
    region = (part.cells,)
    part.SectionAssignment(region=region, sectionName='beam_section', offset=0.0,
                           offsetType=MIDDLE_SURFACE, offsetField='',
                           thicknessAssignment=FROM_SECTION)
    del sketch, region,
    return part


def colum_part(model, colum_width, colum_depth, length, name="colum_part", E=30e9, mu=0.2, material=None, section=None):
    """
    创建柱部件，并且该部件被分配材料与截面
    Parameters
    ----------
    model 模型名称
    colum_width 柱宽度
    colum_depth 柱高度（一般情况与柱宽一致）
    length 柱长度
    name 柱部件名称
    E 弹性模量
    mu 泊松比
    material 材料名称
    section 截面名称

    Returns 创建的柱部件
    -------

    """
    # 创建草图
    sketch = model.ConstrainedSketch(name='colum_sketch', sheetSize=2 * max(colum_width, colum_depth))
    sketch.setPrimaryObject(option=STANDALONE)
    # 绘制草图
    sketch.rectangle(point1=(-colum_width / 2, -colum_depth / 2), point2=(colum_width / 2, colum_depth / 2))
    # 创建部件
    part = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=sketch, depth=length)
    # 如果没有指定该梁的材料或者材料没有被定义，则自定义材料；当材料被定义时，弹性模量E与泊松比mu没用
    if not material:
        material = model.Material(name='concrete')
        material.Elastic(table=((E, mu),))
    # 如果没有指定该梁的截面或者截面为被定义，则自定义截面
    if not section:
        section = model.HomogeneousSolidSection(name='colum_section', material='concrete', thickness=None)
    # 分配截面到实体
    region = (part.cells,)
    part.SectionAssignment(region=region, sectionName='colum_section', offset=0.0,
                           offsetType=MIDDLE_SURFACE, offsetField='',
                           thicknessAssignment=FROM_SECTION)
    del sketch, region,
    return part


def colum_assembly(model, part, loca=(0, 0, 0), name=0):
    """
    为柱装配
    Parameters
    ----------
    model 模型名称
    part 部件名称
    loca 柱定位坐标；初始装配位置为柱的底面中心
    name 装配编号

    Returns 装配后的实体
    -------

    """
    # 创建装配
    name = "colum_part_" + str(name)
    assembly = model.rootAssembly
    instance = assembly.Instance(name=name, part=part, dependent=OFF)
    # 设置装配位置
    assembly.translate(instanceList=(name,), vector=loca)
    return instance


def beam_assembly(model, part, loca=((0, 0), (0, 0)), elevation=3.6, name=0):
    """
    为梁装配
    Parameters
    ----------
    model 模型名称
    part 部件名称
    loca 梁定位坐标；梁需要两个定位坐标：初始端中心点坐标，末端中心点坐标
    elevation 梁顶面高度
    name 装配编号

    Returns 装配后的实体
    -------

    """
    # 创建装配
    name = "beam_part_" + str(name)
    assembly = model.rootAssembly
    instance = assembly.Instance(name=name, part=part, dependent=OFF)
    # 创建参考点与参考平面
    assembly.ReferencePoint(point=(loca[0][0], loca[0][1], elevation))
    ref_point1 = max(assembly.referencePoints.values(), key=lambda x: x)
    assembly.ReferencePoint(point=(loca[1][0], loca[1][1], elevation))
    ref_point2 = max(assembly.referencePoints.values(), key=lambda x: x)
    movable_point1 = instance.InterestingPoint(edge=instance.edges[10], rule=MIDDLE)
    movable_point2 = instance.InterestingPoint(edge=instance.edges[11], rule=MIDDLE)
    assembly.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=3.6)
    ref_plane = max(assembly.datums.values(), key=lambda x: x)
    movable_plane = instance.faces[3]
    # 设定位置
    assembly.CoincidentPoint(movablePoint=movable_point1, fixedPoint=ref_point1)
    assembly.CoincidentPoint(movable_point2, ref_point2)
    assembly.FaceToFace(movable_plane, ref_plane, flip=OFF, clearance=0.0)
    return instance


def set_pressure(model, location, step, magnitude, name=0):
    """
    施加均布荷载
    Parameters
    ----------
    model 模型名称
    location 荷载作用位置；使用findAt搜寻荷载作用的表面
    step 荷载作用的分析步
    magnitude 荷载作用高度
    name 荷载编号

    Returns
    -------

    """
    instance = model.rootAssembly.instances["Structure-1"]
    faces = instance.faces
    for i in range(len(location)):
        side1Faces1 = faces.findAt((location[i],))
        region = model.rootAssembly.Surface(side1Faces=side1Faces1, name="load_" + str(name) + "_" + str(i))
        model.Pressure(name="load_" + str(name) + "_" + str(i), createStepName=step,
                       region=region, distributionType=UNIFORM, field='', magnitude=magnitude,
                       amplitude=UNSET)
    return


def set_restraint(model, location, step, name=0):
    """
    施加固端约束
    Parameters
    ----------
    model 模型名称
    location 固端约束施加的位置
    step 固端约束施加的分析步
    name 固端约束编号

    Returns
    -------

    """
    assembly = model.rootAssembly
    instance = assembly.instances['Structure-1']
    faces = instance.faces
    for i in range(len(location)):
        face = faces.findAt((location[i],))
        region = assembly.Set(faces=face, name="restraint_" + str(name) + "_" + str(i))
        model.DisplacementBC(name="restraint_" + str(name) + "_" + str(i), createStepName=step,
                             region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0,
                             amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='',
                             localCsys=None)
    return


def set_mesh(model, mesh_size=1.0):
    """
    为结构划分网格
    Parameters
    ----------
    model 模型名称
    mesh_size 网格尺寸

    Returns
    -------

    """
    part = model.parts['Structure']
    part.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    part.setMeshControls(regions=part.cells, elemShape=TET, technique=FREE)
    elemType1 = mesh.ElemType(elemCode=C3D20R)
    elemType2 = mesh.ElemType(elemCode=C3D15)
    elemType3 = mesh.ElemType(elemCode=C3D10)
    part.setElementType(regions=(part.cells,), elemTypes=(elemType1, elemType2,
                                                          elemType3))
    part.generateMesh()


def shear_copy(model, number, shear_height=3.6):
    """
    层间复制
    Parameters
    ----------
    model 模型名称
    number 复制的层数
    shear_height 层高

    Returns
    -------

    """
    assembly = model.rootAssembly
    instances = assembly.instances
    instances = tuple([instances.keys()[i] for i in range(len(instances))])
    assembly.LinearInstancePattern(instanceList=instances, direction1=(0.0, 0.0, 1.0), direction2=(0.0, 1.0, 0.0),
                                   number1=number, number2=1, spacing1=shear_height, spacing2=1)
    return


def get_instances(model):
    """
    中间变量处理函数，满足数据传输格式
    """
    assembly = model.rootAssembly
    instances = assembly.instances
    instances_key = tuple([instances.keys()[i] for i in range(len(instances))])
    instances_value = tuple([instances.values()[i] for i in range(len(instances))])
    return instances_key, instances_value


# 创建模型
model = create_model("test")
# 创建梁与柱部件
colum = colum_part(model, 0.4, 0.4, 3.6)
beam = beam_part(model, 0.4, 0.6, 12)
# 创建装配
colums = []
colums.append(colum_assembly(model, colum))
colums.append(colum_assembly(model, colum, (12, 0, 0), 1))
colums.append(colum_assembly(model, colum, (12, 12, 0), 2))
colums.append(colum_assembly(model, colum, (0, 12, 0), 3))
beams = []
beams.append(beam_assembly(model, beam, ((0, 0), (0, 12))))
beams.append(beam_assembly(model, beam, ((0, 12), (12, 12)), 3.6, 1))
beams.append(beam_assembly(model, beam, ((12, 12), (12, 0)), 3.6, 2))
beams.append(beam_assembly(model, beam, ((12, 0), (0, 0)), 3.6, 3))
# 层间复制
shear_copy(model, 8, 3.6)
# 布尔运算为接触部分施加固结约束
instances_key, instances_value = get_instances(model)
model.rootAssembly.InstanceFromBooleanMerge(name='Structure', instances=instances_value, keepIntersections=ON,
                                            originalInstances=SUPPRESS, domain=GEOMETRY)
# 创建分析步
model.StaticStep(name='Step-1', previous='Initial')
# 施加荷载
location_load = []
for i in range(8):
    location_load += [(7.933333, -0.133333, 3.6 * (i + 1)), (12.133333, 7.933333, 3.6 * (i + 1)),
                      (7.933333, 11.866667, 3.6 * (i + 1)), (-0.133333, 4.066667, 3.6 * (i + 1))]
# 施加约束
location_res = (
    (12.066667, -0.066667, 0.0), (0.066667, 11.933333, 0.0), (0.066667, -0.066667, 0.0), (12.066667, 11.933333, 0.0))
set_pressure(model, location_load, "Step-1", 120.0)
set_restraint(model, location_res, "Step-1")
# 划分网格
set_mesh(model, 1.0)
# 创建作业
job = mdb.Job(name='SimpleSupportedBeam_Job', model=model.name, description='Simple Supported Beam Analysis')
job.submit()
mdb.jobs[job.name].waitForCompletion()
# 汇总数据，打开数据库
session.mdbData.summary()
session.openOdb(name='F:/temp/SimpleSupportedBeam_Job.odb')  # 请换成你自己的数据库地址
