# 导入所需模块
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from load import *
from mesh import *
"""
#注释掉的这一部分是ChatGPT给出的结果
# 导入所需模块
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from load import *
from mesh import *

# 创建模型
myModel = mdb.Model(name='CantileverBeamModel')

# 创建实体部件
myPart = myModel.Part(name='CantileverBeamPart', 
                     dimensionality=THREE_D,
                     type=DEFORMABLE_BODY)
mySketch = myPart.ConstrainedSketch(name='BeamProfile', sheetSize=200.0)
mySketch.rectangle(point1=(0.0, 0.0), point2=(0.3, 0.6))
myPart.BaseSolidExtrude(sketch=mySketch, depth=12.0)

# 创建材料
myMaterial = myModel.Material(name='SteelMaterial')
myMaterial.Elastic(table=((209e9, 0.3),))

# 创建截面
mySection = myModel.HomogeneousSolidSection(
name='BeamSection', material='SteelMaterial', thickness=0.6)

# 分配截面属性
myRegion = (myPart.cells,)
myPart.SectionAssignment(region=myRegion, sectionName='BeamSection')

# 装配
myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name='CantileverBeamInstance', 
                               part=myPart, dependent=ON)

# 创建分析步
myStep = myModel.StaticStep(name='LoadStep', previous='Initial')

# 施加均布荷载
pressure = 12000.0  # N/m^2
loadRegion = myPart.faces.findAt(((0.15, 0.3, 6.0),))  # 根据实际情况确定荷载施加的面
myModel.Pressure(name='LoadPressure', 
                createStepName='LoadStep', 
                region=loadRegion, 
                magnitude=pressure)

# 创建网格
elemSize = 1.0  # 网格尺寸
myPart.seedPart(size=elemSize)
myPart.generateMesh()

# 创建作业
myJob = mdb.Job(name='CantileverBeamJob', model='CantileverBeamModel')

# 提交作业
myJob.submit()
myJob.waitForCompletion()

"""
# 创建模型
myModel = mdb.Model(name='CantileverBeamModel')

# 创建实体部件
myPart = myModel.Part(name='CantileverBeamPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)
mySketch = myPart.ConstrainedSketch(name='BeamProfile', sheetSize=200.0)
mySketch.rectangle(point1=(0.0, 0.0), point2=(0.3, 0.6))
myPart.BaseSolidExtrude(sketch=mySketch, depth=12.0)

# 创建材料
myMaterial = myModel.Material(name='SteelMaterial')
myMaterial.Elastic(table=((209e9, 0.3),))

# 创建截面
mySection = myModel.HomogeneousSolidSection(name='BeamSection', material='SteelMaterial', thickness=0.6)

# 分配截面属性
myRegion = (myPart.cells,)
myPart.SectionAssignment(region=myRegion, sectionName='BeamSection')

# 装配
myAssembly = myModel.rootAssembly
myInstance = myAssembly.Instance(name='CantileverBeamInstance', part=myPart, dependent=ON)

# 创建分析步
myStep = myModel.StaticStep(name='LoadStep', previous='Initial')

# 找到荷载施加的面
region_load_center = (0.1, 0.6, 6)
region_load = myInstance.faces.findAt((region_load_center,))
region_load = myAssembly.Surface(side1Faces=region_load, name="load_surface")

# 找到约束施加的面
region_fix_center = (0.1, 0.6, 0)
region_fix = myInstance.faces.findAt((region_fix_center,))
region_fix = myAssembly.Set(faces=region_fix, name="fix_Set")

# 创建均布荷载
myModel.Pressure(name='Load',
                 createStepName='LoadStep', region=region_load, distributionType=UNIFORM,
                 field='', magnitude=12000, amplitude=UNSET)

# 创建固支约束
myModel.EncastreBC(name='FixedSupport', createStepName='LoadStep', region=region_fix)

# 创建网格
elemSize = 1.0  # 网格尺寸
myPart.seedPart(size=elemSize)
myPart.generateMesh()

# 创建作业
myJob = mdb.Job(name='CantileverBeamJob', model='CantileverBeamModel')

# 提交作业
myJob.submit()
myJob.waitForCompletion()
