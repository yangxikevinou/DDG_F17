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
		let A=this.geometry.laplaceMatrix(this.vertexIndex);
		A.scaleBy(h);
		A.incrementBy(M);
		return A;
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		let vertices = this.geometry.mesh.vertices;
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
		// center mesh positions around origin
		normalize(this.geometry.positions, vertices, false);
	}
}
