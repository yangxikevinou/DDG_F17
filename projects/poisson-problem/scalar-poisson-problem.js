"use strict";

/**
 * @module Projects
 */
class ScalarPoissonProblem {
	/**
	 * This class solves a {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf scalar poisson problem} on a surface mesh.
	 * @constructor module:Projects.ScalarPoissonProblem
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @property {number} totalArea The total surface area of the input mesh.
	 */
	constructor(geometry) {
		// index vertices
		this.vertexIndex = indexElements(geometry.mesh.vertices);

<<<<<<< HEAD
<<<<<<< HEAD
		// build the laplace and mass matrices, and compute total area
=======
		// build laplace and mass matrices
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		// build laplace and mass matrices
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		this.A = geometry.laplaceMatrix(this.vertexIndex);
		this.M = geometry.massMatrix(this.vertexIndex);
		this.totalArea = geometry.totalArea();
	}

	/**
	 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A
	 * is the positive definite laplace matrix and M is the mass matrix.
	 * @method module:Projects.ScalarPoissonProblem#solve
	 * @param {module:LinearAlgebra.DenseMatrix} rho A scalar density of vertices of the input mesh.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	solve(rho) {
<<<<<<< HEAD
<<<<<<< HEAD
		let b=this.M.timesDense(rho);
		b=DenseMatrix.constant(b.sum()/this.totalArea,rho.nRows(),1);
		b.decrementBy(rho);
		let llt=this.A.chol();
		return llt.solvePositiveDefinite(this.M.timesDense(b));
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		// construct right hand side
		let V = this.M.nRows();
		let totalRho = this.M.timesDense(rho).sum();
		let rhoBar = DenseMatrix.ones(V, 1).timesReal(totalRho / this.totalArea);
		let rhs = this.M.timesDense(rhoBar.minus(rho));

		// solve linear system
		let llt = this.A.chol();
		let phi = llt.solvePositiveDefinite(rhs);

		return phi;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
	}
}
