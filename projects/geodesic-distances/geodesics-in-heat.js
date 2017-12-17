"use strict";

class GeodesicsInHeat {
	/**
	 * This class implements the {@link https://www.cs.cmu.edu/~kmcrane/Projects/GeodesicsInHeat geodesics in heat} algorithm to compute geodesic distances
	 * on a surface mesh.
	 * @constructor module:Projects.GeodesicsInHeat
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} F The mean curvature flow operator built on the input mesh.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);

		this.A = geometry.laplaceMatrix(this.vertexIndex);
		this.F = geometry.massMatrix(this.vertexIndex).plus(this.A.timesReal(Math.pow(geometry.meanEdgeLength(),2)));
	}

	/**
	 * Computes the vector field X = -∇u / |∇u|.
	 * @private
	 * @method module:Projects.GeodesicsInHeat#computeVectorField
	 * @param {module:LinearAlgebra.DenseMatrix} u A dense vector (i.e., u.nCols() == 1) representing the
	 * heat that is allowed to diffuse on the input mesh for a brief period of time.
	 * @returns {Object} A dictionary mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 */
	computeVectorField(u) {
		let X={};
		let G=this.geometry;
		for (let f of G.mesh.faces) {
			if (!f.isBoundaryLoop()) {
				let x=new Vector();
				for (let h of f.adjacentHalfedges()) {
					x.decrementBy(G.vector(h).times(u.get(this.vertexIndex[h.prev.vertex],0)));
				}
				x=G.faceNormal(f).cross(x);
				x.normalize();
				X[f]=x;
			}
		}
		return X;
	}

	/**
	 * Computes the integrated divergence ∇.X.
	 * @private
	 * @method module:Projects.GeodesicsInHeat#computeDivergence
	 * @param {Object} X The vector field -∇u / |∇u| represented by a dictionary
	 * mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeDivergence(X) {		
		let G=this.geometry;
		let VS=G.mesh.vertices;
		let div=DenseMatrix.zeros(VS.length,1);
		for (let v of VS) {
			let D=0.;
			for (let h of v.adjacentHalfedges()) {
				D+=G.cotan(h)*(G.vector(h).dot(X[h.face]));
				h=h.twin;
				D-=G.cotan(h)*(G.vector(h).dot(X[h.face]));
			}
			div.set(D/2,this.vertexIndex[v],0);
		}
		return div;
	}

	/**
	 * Shifts φ such that its minimum value is zero.
	 * @private
	 * @method module:Projects.GeodesicsInHeat#subtractMinimumDistance
	 * @param {module:LinearAlgebra.DenseMatrix} phi The (minimum 0) solution to the poisson equation Δφ = ∇.X.
	 */
	subtractMinimumDistance(phi) {
		let min = Infinity;
		for (let i = 0; i < phi.nRows(); i++) {
			min = Math.min(phi.get(i, 0), min);
		}

		for (let i = 0; i < phi.nRows(); i++) {
			phi.set(phi.get(i, 0) - min, i, 0);
		}
	}

	/**
	 * Computes the geodesic distances φ using the heat method.
	 * @method module:Projects.GeodesicsInHeat#compute
	 * @param {module:LinearAlgebra.DenseMatrix} delta A dense vector (i.e., delta.nCols() == 1) containing
	 * heat sources, i.e., u0 = δ(x).
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	compute(delta) {
		let Ffac=this.F.chol();
		let Afac=this.A.chol();
		let phi=delta;
		phi=this.computeVectorField(Ffac.solvePositiveDefinite(phi));
		phi=Afac.solvePositiveDefinite(this.computeDivergence(phi).timesReal(-1.));		

		// since φ is unique up to an additive constant, it should
		// be shifted such that the smallest distance is zero
		this.subtractMinimumDistance(phi);

		return phi;
	}
}