from bifenics import BifenicsProblem, ParameterContinuation, evaluate_function
from dolfin import (
    Constant,
    Mesh,
    FunctionSpace,
    VectorElement,
    FiniteElement,
    MixedElement,
    SubDomain,
    near,
    Function,
    MeshFunction,
    DirichletBC,
    Measure,
    derivative,
    grad,
    sqrt,
    dot,
    split,
    FacetNormal,
    dx,
    tr,
    det,
    SpatialCoordinate,
    CompiledSubDomain,
    Identity,
    project,
)
from ufl import cofac, sign

import scipy.io
import numpy as np


class DiscoTensioneSuperficiale(BifenicsProblem):

    # Subdomains
    class Symmetry(SubDomain):
        def __init__(self, depth):
            self.depth = depth
            SubDomain.__init__(self)

        def inside(self, x, on_boundary):
            TOL = 1e-10
            return on_boundary and near(x[0], 0, TOL) and x[1] < 2 - self.depth + 1e-5

    def __init__(self, mu=1, depth=1):
        # Elastic constants
        self.mu = Constant(mu)
        self.depth = depth
        self.export_data = {
            "gamma": np.array([]),
            "opening": np.array([]),
        }

    def mesh(self):
        mesh = Mesh("mesh.xml")
        return mesh

    def function_space(self, mesh):
        P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
        P2 = VectorElement("CG", mesh.ufl_cell(), 2)
        THelem = MixedElement(P2, P1)
        V = FunctionSpace(mesh, THelem)
        return V

    def parameters(self):
        return {"gamma": Constant(0)}

    def residual(self, up, vq, parameters):
        gamma = parameters["gamma"]
        mesh = up.function_space().mesh()
        u, p = split(up)
        X = SpatialCoordinate(mesh)
        R = sqrt(X[0] * X[0] + (X[1] - 1) * (X[1] - 1))
        boundaries = self.boundary_function(mesh)
        my_ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

        X = SpatialCoordinate(mesh)
        x = X + u

        F = grad(x)
        J = det(F)
        C = F.T * F
        I1 = tr(C) * J ** (-2.0 / 3.0)
        Wbulk = (0.5 * self.mu * (I1 - 2) - p * (J - 1)) * dx

        N = FacetNormal(mesh)
        NansonOp = cofac(F)
        deformed_N = dot(NansonOp, N)
        current_element_of_area = sqrt(dot(deformed_N, deformed_N))

        surface_energy_density = gamma * current_element_of_area
        Wsurface = surface_energy_density * my_ds(0)

        W = Wbulk + Wsurface

        FF = derivative(W, up, vq)
        return FF

    def solver_parameters(self):
        solver_params = {
            "nonlinear_solver": "snes",
            "snes_solver": {
                "linear_solver": "mumps",
                "absolute_tolerance": 1e-8,
                "relative_tolerance": 1e-8,
                "maximum_iterations": 10,
                "error_on_nonconvergence": False,
            },
        }

        return solver_params

    def boundary_function(self, mesh):
        symmetry = self.Symmetry(self.depth)

        boundaries = MeshFunction("size_t", mesh, 1)
        boundaries.set_all(0)

        symmetry.mark(boundaries, 1)
        return boundaries

    def boundary_conditions(self, mesh, V):
        boundaries = self.boundary_function(mesh)

        corner = CompiledSubDomain("near(x[1], 0, TOL) && near(x[0], 0, TOL)", TOL=1e-6)

        bcc = DirichletBC(V.sub(0), Constant((0, 0)), corner, method="pointwise")
        bcr = DirichletBC(V.sub(0).sub(0), Constant(0), boundaries, 1)

        return [bcc, bcr]

    def monitor(self, up, parameters, xdmf_file):
        gamma = parameters["gamma"]
        gamma_float = round(float(gamma), 10)

        mesh = up.function_space().mesh()
        comm = mesh.mpi_comm()
        u = Function(up, 0, name="displacement")
        p = Function(up, 1, name="pressure")

        t = float(gamma)
        t = round(abs(t), 10)

        with xdmf_file as xdmf:
            xdmf.write(u, t)
        with xdmf_file as xdmf:
            xdmf.write(p, t)

        spos = evaluate_function(u, (0, 2))[0]

        self.export_data["gamma"] = np.append(self.export_data["gamma"], gamma_float)
        self.export_data["opening"] = np.append(self.export_data["opening"], spos)

        scipy.io.savemat("output/data.mat", self.export_data)

        Id = Identity(2)
        F = Id + grad(u)
        B = F * F.T

        tr_cauchy = tr(self.mu * (B - Id) - p * Id)

        W = FunctionSpace(mesh, "CG", 1)
        interpolated_cauchy = project(tr_cauchy, W)
        interpolated_cauchy.rename("cauchy", "cauchy")

        xdmf_file.write(interpolated_cauchy, t)


if __name__ == "__main__":
    disco_tens_sup = DiscoTensioneSuperficiale(mu=1, depth=1.2)
    XDMF_options = {
        "flush_output": True,
        "functions_share_mesh": True,
        "rewrite_function_mesh": False,
    }
    analysis = ParameterContinuation(
        disco_tens_sup,
        "gamma",
        start=0,
        end=0.2,
        dt=0.002,
        save_output=False,
        saving_file_parameters=XDMF_options,
    )
    analysis.run()
