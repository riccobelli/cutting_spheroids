import gmsh
import os

gmsh.initialize()

gmsh.model.add("tumor_cut")

ri = 0.04
ro = 1
depth = 1

gmsh.model.occ.addSphere(0, 0, 0, ro, 1)
M = 1
cut = []
gmsh.model.occ.addBox(-ro, -ro, -ro, 2 * ro, ro, 2 * ro, M + 1)
gmsh.model.occ.addBox(-ro, -ro, -ro, ro, 2 * ro, 2 * ro, M + 2)
gmsh.model.occ.addCylinder(-ro, 0, (1 - depth + ri), 2 * ro, 0, 0, ri, M + 4)
gmsh.model.occ.addBox(0, 0, (1 - depth + ri), ro, ri, 2 * ro, M + 3)
for i in range(1, 5):
    cut.append((3, M + i))
ov, _ = gmsh.model.occ.cut([(3, 1)], cut)
gmsh.model.occ.synchronize()

gmsh.model.mesh.field.add("MathEval", 1)
gmsh.model.mesh.field.setString(
    1, "F", f"0.1*(1 - exp(-100*(y^2+(z-(1-{depth}+{ri}))^2)))"
)
gmsh.model.mesh.field.add("MathEval", 2)
gmsh.model.mesh.field.setString(
    2,
    "F",
    f"0.1*(1 - exp(-50*(x^2+z^2-1)^2-100*y^2))*(tanh(100*(z-(1-{depth}+{ri})))+1)/2+0.1*(tanh(-100*(z-(1-{depth}+{ri})))+1)/2",
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

gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("mesh.msh")

gmsh.finalize()
os.system("dolfin-convert mesh.msh mesh.xml")
