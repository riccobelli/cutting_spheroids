import gmsh
import os

gmsh.initialize()

gmsh.model.add("tumor_cut")

ro = 1
depth = 1.2

gmsh.model.occ.addDisk(0, 1, 0, 1, 1, 1)
M = 1
cut = []
gmsh.model.occ.addRectangle(-1, 0, 0, 1, 2, 2)
ov, _ = gmsh.model.occ.cut([(2, 1)], [(2, 2)])
gmsh.model.occ.synchronize()

gmsh.model.mesh.field.add("MathEval", 1)
gmsh.model.mesh.field.setString(1, "F", f"0.05*(1 - exp(-100*(x^2+(y-2+{depth})^2)))")
gmsh.model.mesh.field.add("MathEval", 2)
gmsh.model.mesh.field.setString(
    2,
    "F",
    "0.05*(1 - exp(-100*(x^2+(y-2)^2)))",
)
gmsh.model.mesh.field.add("MathEval", 3)
gmsh.model.mesh.field.setString(3, "F", "0.005")

gmsh.model.mesh.field.add("Max", 4)
gmsh.model.mesh.field.setNumbers(4, "FieldsList", [1, 3])
gmsh.model.mesh.field.add("Max", 5)
gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 3])
gmsh.model.mesh.field.add("Min", 6)
gmsh.model.mesh.field.setNumbers(6, "FieldsList", [4, 5])

gmsh.model.mesh.field.setAsBackgroundMesh(6)


gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.option.setNumber("Mesh.Algorithm", 5)

gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("mesh.msh")

gmsh.finalize()
os.system("dolfin-convert mesh.msh mesh.xml")
