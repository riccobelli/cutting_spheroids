from bifenics import BifenicsProblem, ParameterContinuation, evaluate_function
from dolfin import (
    Constant,
    Mesh,
    SubDomain,
    near,
    MeshFunction,
    DirichletBC,
    Measure,
    derivative,
    grad,
    sqrt,
    sin,
    cos,
    conditional,
    dot,
    FacetNormal,
    dx,
    tr,
    det,
    SpatialCoordinate,
    CompiledSubDomain,
    Identity,
    FunctionSpace,
    parameters,
    as_vector,
    outer,
    FiniteElement,
    VectorElement,
    MixedElement,
    split,
    Function,
)
from ufl import cofac, atan_2
import os
import scipy.io
import numpy as np

parameters["form_compiler"]["quadrature_degree"] = 6
os.environ["OMP_NUM_THREADS"] = "1"


class SferaTensioneSuperficiale(BifenicsProblem):

    # Subdomains
    class xz_Non_cut(SubDomain):
        def inside(self, x, on_boundary):
            TOL = 1e-10
            return on_boundary and near(x[1], 0, TOL)

    class yz_Plane(SubDomain):
        def inside(self, x, on_boundary):
            TOL = 1e-10
            return on_boundary and near(x[0], 0, TOL)

    def __init__(self, mu=1, simulate_growth=False):
        self.mu = Constant(mu)
        self.simulate_growth = simulate_growth
        self.export_data = {
            "gamma": np.array([]),
            "opening": np.array([]),
        }

    def mesh(self):
        mesh = Mesh("mesh.xml")
        return mesh

    def function_space(self, mesh):
        VP2elem = VectorElement("CG", mesh.ufl_cell(), 2)
        P0elem = FiniteElement("DG", mesh.ufl_cell(), 0)
        elem = MixedElement([VP2elem, P0elem])
        V = FunctionSpace(mesh, elem)
        return V

    def parameters(self):
        return {"gamma": Constant(0)}

    def residual(self, up, vq, parameters):
        gamma = parameters["gamma"]
        mesh = up.function_space().mesh()
        u, p = split(up)

        boundaries = self.boundary_function(mesh)
        my_ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
        if self.simulate_growth == True:
            X = SpatialCoordinate(mesh)
            R = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2])
            PHI = atan_2(X[1], X[0])
            RHO = sqrt(X[0] ** 2 + X[1] ** 2)
            THETA = atan_2(RHO, X[2])
            sinT, cosT = [sin(THETA), cos(THETA)]
            sinP, cosP = [sin(PHI), cos(PHI)]
            ER = as_vector([sinT * cosP, sinT * sinP, cosT])
            ET = as_vector([cosT * cosP, cosT * sinP, -sinT])
            EP = as_vector([-sinP, cosP, 0])
            tmp = 96 * R**3 * gamma**3 * self.mu**3 + 81 * R**6 * self.mu**6
            tmp2 = 9 * R**3 * self.mu**3 + sqrt(tmp)
            pR = (
                -4 * 3 ** (1 / 3) * R * gamma * self.mu
                + 2 ** (1 / 3) * (tmp2) ** (2 / 3)
            ) / (6 ** (2 / 3) * R * self.mu * (tmp2) ** (1 / 3))
            Pgrowth = pR * outer(ER, ER) + 1 / sqrt(pR) * (Identity(3) - outer(ER, ER))
            P = conditional(R > 0.5, Pgrowth, Identity(3))
        else:
            P = Identity(3)
        X = SpatialCoordinate(mesh)
        x = X + u

        F = grad(x)
        J = det(F)
        Fe = F * P
        Je = det(Fe)
        Ce = Fe.T * Fe
        I1e = tr(Ce)
        I1ebar = I1e / Je ** (2 / 3)
        Wbulk = (0.5 * self.mu * (I1ebar - 3) - p * (J - 1)) * dx

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
                "absolute_tolerance": 1e-10,
                "relative_tolerance": 1e-10,
                "maximum_iterations": 10,
                "error_on_nonconvergence": False,
            },
        }

        return solver_params

    def boundary_function(self, mesh):
        xz_non_cut = self.xz_Non_cut()
        yz_plane = self.yz_Plane()

        boundaries = MeshFunction("size_t", mesh, 2)
        boundaries.set_all(0)

        xz_non_cut.mark(boundaries, 1)
        yz_plane.mark(boundaries, 2)
        return boundaries

    def boundary_conditions(self, mesh, V):
        boundaries = self.boundary_function(mesh)

        corner = CompiledSubDomain(
            "near(x[0], 0, TOL) && near(x[1], 0, TOL) && near(x[2], -1, TOL)", TOL=1e-6
        )

        bcc = DirichletBC(V.sub(0), Constant((0, 0, 0)), corner, method="pointwise")
        bcxz = DirichletBC(V.sub(0).sub(1), Constant(0), boundaries, 1)
        bcyz = DirichletBC(V.sub(0).sub(0), Constant(0), boundaries, 2)

        return [bcc, bcxz, bcyz]

    def monitor(self, up, parameters, xdmf_file):
        gamma = parameters["gamma"]
        gamma_float = round(float(gamma), 10)

        mesh = up.function_space().mesh()
        u = Function(up, 0, name="displacement")
        p = Function(up, 1, name="pressure")
        point = max(mesh.coordinates()[:], key=lambda x: x[2])
        comm = mesh.mpi_comm()
        points = comm.gather(point, root=0)
        if comm.rank == 0:
            point = max(points, key=lambda x: x[2])
        point = comm.bcast(point, root=0)
        spos = evaluate_function(u, point)[1]

        self.export_data["gamma"] = np.append(self.export_data["gamma"], gamma_float)
        self.export_data["opening"] = np.append(self.export_data["opening"], spos)
        scipy.io.savemat("output/data.mat", self.export_data)

        t = float(gamma)
        t = round(abs(t), 10)
        with xdmf_file as xdmf:
            xdmf.write(u, t)
        with xdmf_file as xdmf:
            xdmf.write(p, t)


if __name__ == "__main__":
    sfera_tens_sup = SferaTensioneSuperficiale(mu=1, simulate_growth=False)
    XDMF_options = {
        "flush_output": True,
        "functions_share_mesh": True,
        "rewrite_function_mesh": False,
    }
    analysis = ParameterContinuation(
        sfera_tens_sup,
        "gamma",
        start=0,
        end=0.2,
        dt=0.01,
        saving_file_parameters=XDMF_options,
        max_dt=0.01,
        max_step_for_dt_doubling=3,
        save_output=False,
    )
    analysis.run()
