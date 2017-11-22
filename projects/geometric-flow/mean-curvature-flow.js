"use strict";

class MeanCurvatureFlow {
	/**
	 * This class performs {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the mean curvature flow operator.
	 * @private
	 * @method module:Projects.MeanCurvatureFlow#buildFlowOperator
	 * @param {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @param {number} h The timestep.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	buildFlowOperator(M, h) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let A=this.geometry.laplaceMatrix(this.vertexIndex);
		A.scaleBy(h);
		A.incrementBy(M);
		return A;
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let A = this.geometry.laplaceMatrix(this.vertexIndex);

		// F = M + hA
		return M.plus(A.timesReal(h));
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		// build the flow and mass matrices
		let vertices = this.geometry.mesh.vertices;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let f=DenseMatrix.zeros(vertices.length,3);
		for (let v of vertices) {
			let p=this.geometry.positions[v];
			f.set(p.x,this.vertexIndex[v],0);
			f.set(p.y,this.vertexIndex[v],1);
			f.set(p.z,this.vertexIndex[v],2);
		}
		let M=this.geometry.massMatrix(this.vertexIndex);
		let llt=this.buildFlowOperator(M,h).chol();
		let g=llt.solvePositiveDefinite(M.timesDense(f));
		for (let v of vertices) {
			this.geometry.positions[v]=new Vector(g.get(this.vertexIndex[v],0),g.get(this.vertexIndex[v],1),g.get(this.vertexIndex[v],2));
		}
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let V = vertices.length;
		let M = this.geometry.massMatrix(this.vertexIndex);
		let F = this.buildFlowOperator(M, h);

		// construct right hand side
		let f0 = DenseMatrix.zeros(V, 3);
		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let p = this.geometry.positions[v];
<<<<<<< HEAD
=======

			f0.set(p.x, i, 0);
			f0.set(p.y, i, 1);
			f0.set(p.z, i, 2);
		}

		let rhs = M.timesDense(f0);

		// solve linear system (M - hA)fh = Mf0
		let llt = F.chol();
		let fh = llt.solvePositiveDefinite(rhs);

		// update positions
		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let p = this.geometry.positions[v];

			p.x = fh.get(i, 0);
			p.y = fh.get(i, 1);
			p.z = fh.get(i, 2);
		}
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5

			f0.set(p.x, i, 0);
			f0.set(p.y, i, 1);
			f0.set(p.z, i, 2);
		}

		let rhs = M.timesDense(f0);

		// solve linear system (M - hA)fh = Mf0
		let llt = F.chol();
		let fh = llt.solvePositiveDefinite(rhs);

		// update positions
		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let p = this.geometry.positions[v];

			p.x = fh.get(i, 0);
			p.y = fh.get(i, 1);
			p.z = fh.get(i, 2);
		}

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		// center mesh positions around origin
		normalize(this.geometry.positions, vertices, false);
	}
}
